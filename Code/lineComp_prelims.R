# SETUP ----------------------------
# ~ Code 
library(multitaper) # spec.mtm, etc.
library(pbapply)    # Progress bar for apply() functions
library(itsmr)      # Time series stuff, used here for ARMA
library(nleqslv)    # Stupidly named nonlinear solving package
library(nloptr)     # A more flexible nonlinear solver
library(JADE)       # BSS

# ~ Presentation (Plots and  Tables) 
library(kableExtra) # Nice Tables
library(animation)  # For creating gifs
library(fields)     # Supplement image plots: legends, better colour schema

load("~/Skye_Toolbox/Splot/splot.RData")
load("~/Skye_Toolbox/Splot/spalette.RData")

# ~ Options 
options(digits = 10)
# --------------------------------------


# PRELIMS ------------------------------
# ~ Constants & Multitaper Parameters ----
set.seed(68)
i  <- complex(1,0,1)          # imaginary i

N  <- 1000                    # No. of Observations
NF <- 2^ceiling(log(N,2))+1   # No. of Fourier freqs
FF <- 1/(2*NF)                # Fundamental Fourier freq
f  <- FF*(0:(NF-1))           # Freq vector
t  <- 1:N                     # Time vector

w  <- 4/N                     # Analysis bandwidth (freq)
K  <- 7                       # Eigencoefs per frequency
v  <- spec.mtm(ts(t),plot = F)$mtm$dpss$v    # DPSS(N)

# ~ Line Components ----
L <- 5
freqs   <- seq(0.025,0.475, length.out = L)
amps    <- 1:L*0.25
signals <- matrix(0, N, L)
for(j in 1:L){signals[,j] <- amps[j]*cos(2*pi*freqs[j]*1:1000)}
s.t <- ts(rowSums(signals))

# Quick function for plotting line comp location
sigs <- function(col = "green3"){abline(v = freqs, col = col, lty = 3)}

# ~ c(t): modulating functions ----
ct <- list()
ct$const  <- rep(1,N)

ct$grow   <- list(line  = 0.5 + t/N,
                  sqrt  = ((0.1 + t/100)^(1/2))/2,
                  cube  = 0.5 + (t/600)*((t/600)-1)^2 
)
ct$decay  <- lapply(ct$grow, rev)

ct$center <- list(norm  = exp((-(t-N/2)^2)/(2*(N/5)^2)),
                  boost = exp((-(t-N/2)^2)/(2*(N/5)^2)) + 0.5
)
ct$bounds <- list(norm  = 2 - ct$center$boost,
                  boost = 2 - ct$center$norm
)
ct$wiggle <- 0.5 + (t/800)*((t/800)-1)^2

# g(t): squared modulating functions
gt <- lapply(ct, function(x)lapply(x, function(y)y^2))
gt$const  <- unlist(gt$const)
gt$wiggle <- unlist(gt$wiggle)

# ~ Y(t): stationary time series ----
Y.funk <- function(phi = 0.5, SD = 1, N = 1000){
  ts(arima.sim(model = list(ar = phi, sd = SD), n = N))
}

# ~ Sy(f): stationary spectrum (AR only) ----
AR.spec <- function(phi, sd = 1){
  denom <- 1
  for(h in 1:length(phi)){
    denom <- denom - phi[h]*exp(-i*2*pi*f*h)  # this is 1 - sum(e^i{tau}fh)
  }
  return( sd^2/Mod(denom)^2 )
}

# --------------------------------------


# MATH ---------------------------------
# Standardize (Glen mode) ----
normy <- function(x){ (x-mean(x))/sqrt(var(x)) }

# ~ Create UMP and store its components ----
create <- function(base = ts(rnorm(N)), c.t = ct$const, s.t = s.t, sigs.mod = TRUE){
  
  # Build, considering optional modulation of line components
  if(sigs.mod){ y.t <- base + s.t } else { y.t <- base + s.t/c.t }
  
  # Store and return
  UMP <- list(c = ts(c.t), y = ts(y.t), x = ts(c.t*y.t))
  # return(UMP)
  return(ts(c.t*y.t))
}

# ~ Create pbc ----
pbc  <- function(phi = phii, ct.1 = ct.0, ct.2 = 0, s.t = 0, sigs.mod = TRUE,
                 B = 100, M = 1, sd = 1, seed = 0, p = 1:2){
  
  # B-dependent Prelims 
  B2 <- ceiling(B/2)                                  # Midpoint shift
  b  <- 1:(N-(B-1))                                   # block indices
  tb <- b + B2 - 1                                    # time in-block
  vb <- t(spec.mtm(ts(1:B), plot = FALSE)$mtm$dpss$v) # DPSS (based on B)
  BS <- list(B = B, B2 = B2, b = b, tb = tb, vb = vb) # store for output
  
  
  # ~ Math ----------------
  # N x inner sum coef                    /----w_B----\
  Ma  <- function(t,m){return(B*(sin(2*pi *(K+1)/(2*B)* (t-m))/(pi*(t-m)))^2)}
  
  # Matrix of Interest
  mat <- matrix(0,nrow = B, ncol = B)
  for(s in 1:B){mat[s,] <- Ma(s,1:B)}
  
  # Fill diagonal with Sinc function at "0/0"
  diag(mat) <- 1 
  
  # Get eigenstuff
  mval <- eigen(mat)$values
  mvec <- eigen(mat)$vectors
  
  # Function to estimate coefficients {a(l,b)}
  a  <- function(s,l,b){ (K/(B*(mval[l+1]))) * (matrix(s[,b], nrow = NF) %*% mvec[,l+1]) }
  
  # omega function vals for upcoming Taylor expansions (before loop for speed)
  w0.0 <- sum(mvec[,1])
  w2.0 <- sum(mvec[,3])
  w1.1 <- c(((-B2):(B2 - 1 - (B/2 != B2)))   %*% mvec[,2])
  w0.2 <- c(((-B2):(B2 - 1 - (B/2 != B2)))^2 %*% mvec[,1])
  w2.2 <- c(((-B2):(B2 - 1 - (B/2 != B2)))^2 %*% mvec[,3])
  d.02 <- (w2.0*w0.2 - w0.0*w2.2)
  
  # Functions for first and second time-derivatives
  tds  <- function(s,p){ B*a(s,1,p)/w1.1 }
  stds <- function(s,p){ 2*B * (a(s,0,p)*w2.0 - a(s,2,p)*w0.0) / d.02}
  
  
  # ~ LOOP ----
  
  # Loop setup 
  if(seed > 0){ set.seed(seed) }
  pb <- txtProgressBar(style = 3)
  
  # Initialize arrays
  bc1   <- if(any(p==1)) {array(dim = c(NF,N,M))} else { NULL }
  bc2   <- if(any(p==2)) {array(dim = c(NF,N,M))} else { NULL }
  ##ftest <- array(0, dim = c(NF,N,M))
  
  # ...........................................................
  for(m in 1:M){
    
    xt1 <- create(base = Y.funk(phi), c.t = ct.1, s.t = s.t, sigs.mod = sigs.mod)
    xt2 <- if(ct.2 != 0){
      create(base = rnorm(N), c.t = ct.2, s.t = s.t, sigs.mod = sigs.mod)
      } else {0}
    xt <- xt1 + xt2
    
    b.list <- list()
    for(bb in b){b.list[bb] <- bb}
    
    # ~ Math ----
    ## sw.mtm <- lapply(b.list, function(bb){
    ##       spec.mtm(ts(xt[(bb):(bb+B-1)]),
    ##                nFFT = (2*NF-1),
    ##                plot = FALSE,
    ##                returnInternals = TRUE,
    ##                Ftest = TRUE)})
    ## 
    # Sliding window: full window included
    ## sw <- (B/K) * Mod( sapply(b, function(bb){ sw.mtm[[bb]]$mtm$eigenCoefs %*% vb }))^2 
    
    sw <- (B/K) * Mod( sapply(b, function(bb){
      (spec.mtm(ts(xt[(bb):(bb+B-1)]),
                nFFT = (2*NF-1),
                plot = FALSE,
                returnInternals = TRUE)$mtm$eigenCoefs) %*% vb }))^2 
    
    # get Ftests
    ##ftest[,tb,m] <- sapply(b, function(bb){sw.mtm[[bb]]$mtm$Ftest} )
    
    # Endpoints: 1st derivative
    tds.start  <- tds(sw,1)
    tds.end    <- tds(sw,max(b)) 
    
    # Endpoints: 2nd derivative
    stds.start <- stds(sw,1)      
    stds.end   <- stds(sw,max(b))
    
    # Heavy lifting done 
    gc(); setTxtProgressBar(pb, (m-0.5)/M)
    
    
    # ~ 1-BC ----
    if(any(p==1)){
      
      # Modify (linear) storing midpoints
      bc1[,tb,m] <- (K/(mval[1]*w0.0)) * 
        sapply(b, function(p){matrix(sw[,p], nrow = NF) %*% mvec[,1]})
      
      # Extrapolate
      bc1[,  1:(tb[1]-1),m] <- sapply((tb[1]-1):1,
                                      function(h){ bc1[,tb[1],m] - h*tds.start })
      bc1[,(max(tb)+1):N,m] <- sapply(1:(B2-(B2 != B/2)),
                                      function(h){ bc1[,max(tb),m] + h*tds.end })
    }
    
    
    # ~ 2-BC ----
    if(any(p==2)){
      
      # Modify (quadratic) storing midpoints
      bc2[,tb,m] <- (B/d.02) * sapply(b, function(p){ (a(sw,2,p)*w0.2 - a(sw,0,p)*w2.2) })
      
      # Extrapolate
      bc2[,  1:(tb[1]-1),m] <- sapply((tb[1]-1):1, 
                                      function(h){ 
                                        bc2[,tb[1],m] - h*tds.start - (h^2)*stds.start/2})
      bc2[,(max(tb)+1):N,m] <- sapply(1:(B2-(B2 != B/2)), 
                                      function(h){ 
                                        bc2[,max(tb),m] + h*tds.end + (h^2)*stds.end/2  })
    }
    
    # Clean
    rm(sw); gc()
    
    # Phew!
    setTxtProgressBar(pb, m/M); gc()
  } ## --------------------------- loop ends --------------------------- ##
  
  # output
  return(list(bc1 = bc1,
              bc2 = bc2,
              ## Ftest = ftest,
              B.vars = BS,
              xt.last = list(xt = xt, xt1 = xt1, xt2 = xt2)))
}

# ~ Mean PBC function ----
pbc.means <- function(PBC){
  return(list(
    bc1   = pbapply(PBC$bc1,   1:2, mean)##,
    # bc2   = pbapply(PBC$bc2,   1:2, mean)##,
    ##Ftest = pbapply(PBC$Ftest,   1:2, mean)
  ) )
}

# ~ Smooth g(t) recovery ----
gest <- function(sgram){
  N <- length(sgram[,1])
  A <- outer(1:N,1:N, function(l,j){rowSums(sgram)[l]/rowSums(sgram)[j]})
  G <- Mod(eigen(A/N)$vectors[,1])
  return(G); gc()
}

# ~ Smooth Sy(f) recovery ----
yest <- function(sgram){
  N <- length(sgram[1,])
  A <- outer(1:NF,1:NF, function(l,j){colSums(sgram)[l]/colSums(sgram)[j]})
  S <- Mod(eigen(A/NF)$vectors[,1])
  return(S); gc()
}

# ~ SVD recoveries ----
est.a <- function(sgram){
  # SVD
  duv <- svd(sgram,1,1); gc()
  
  d <- duv$d
  u <- duv$u
  v <- duv$v
  
  # Normalize
  bu <- 1/duv$u[1]
  bv <- duv$d[1]/bu
  
  # Estimate g
  g <- bu * duv$u
  
  # Estimate Sy
  sy <- bv * duv$v
  
  # Return
  return(list(g = c(g), s = c(sy)))
}

# ~ Decompose spectrogram ----
decompose <- function(sgram, normalize = TRUE,
                      gg = TRUE, ss = TRUE, gg.svd = TRUE, ss.svd = TRUE){
  
  g     <- if(gg)    {if(normalize){ normy( gest(sgram)  )} else{  gest(sgram)  }}
  s     <- if(ss)    {if(normalize){ normy( yest(sgram)  )} else{  yest(sgram)  }}
  svd.g <- if(gg.svd){if(normalize){ normy(est.a(sgram)$g)} else{ est.a(sgram)$g}}
  svd.s <- if(ss.svd){if(normalize){ normy(est.a(sgram)$s)} else{ est.a(sgram)$s}}
  
  return(list(g = g, s = s, svd.g = svd.g, svd.s = svd.s))
}

# ~ Create smooth TFS estimate, given component vecs ----
recompose <- function(G1, s1, G2 = NULL, s2 = NULL){ 
  if(is.null(G2) || is.null(s2)){ 
    # UMP
    return(outer(G1,s1))
  } else {
    # GUMP-2
    return(outer(G1,s1) + outer(G2,s2))  
  }
}


# ~ Partition-column Smoothing ----
gester <- function(sgram, n1, n2){
  N <- length(sgram[,1])
  A <- outer(1:N, 1:N, function(l,j){
    rowSums(sgram[,n1:n2])[l]/rowSums(sgram[,n1:n2])[j]
  })
  G <- Mod(eigen(A/N)$vectors[,1])
  return(G); gc()
}

jester <- function(sgram, c.t = 1, n = 50, plot = TRUE, mod = 1){
  pb <- txtProgressBar(style = 3)
  
  N   <- length(sgram[,1])
  fit <- N %/% n
  out <- matrix(0, N, fit)
  svd.g <- normy(est.a(sgram)$g)
  old.g <- normy(gest(sgram))
  
  if(plot){
    
    plot.true(c.t^2, mod = mod)
    for(m in 1:fit){
      out[,m] <- gester(sgram,((m-1)*n)+1, m*n)
      lines(t, normy(out[,m]), col = "#0000FF33")  
      setTxtProgressBar(pb, m/fit)
    } 
    if(N%%n != 0){ out <- cbind(out, gester(sgram,m*n+1,N)) }
    
    lines(t, old.g, col = "goldenrod", lwd = 2)
    lines(t, svd.g, col = "red", lwd = 2)
    lines(t, normy(rowMeans(out)), col = "dodgerblue", lwd = 2)
    legend("topright", horiz = TRUE,
           c("true","ours (new)","E[ours (new)]", "ours (old)", "SVD"),
           col = c("black","blue","dodgerblue","goldenrod","red"),
           lty = 1, cex = 0.8, lwd = c(1,1,2,2,2))
    
    # ERROR
    MSE.n <- mean(abs(normy(rowMeans(out)) - normy(c.t^2)))
    mtext(paste("n =",n,"; MSE =", MSE.n))
    
  } else {
    for(m in 1:fit){
      out[,m] <- gester(sgram,((m-1)*n)+1, m*n)
      setTxtProgressBar(pb, m/fit)
    }
    if(N%%n != 0){ out <- cbind(out, gester(sgram,m*n+1,N)) }
    
    # ERROR
    MSE.n <- mean(abs(normy(rowMeans(out)) - normy(c.t^2)))
  }
  
  return(list(out = out, mean = normy(rowMeans(out)), mse = MSE.n,
              old = old.g, svd = svd.g))
}