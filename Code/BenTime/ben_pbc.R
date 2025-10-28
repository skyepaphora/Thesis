p.bc <- function(p = 1, type = "noise", B = 100, M = 1, seed = 1, sd = 10,
                 ct1 = rep(1,N),      ct2 = rep(1,N), 
                 ar2 = c(0.1,-0.7), ar2.2 = c(0.5,-0.5),
                 freq = NULL){
  
  
  # B-dependent Prelims
  B2 <- ceiling(B/2)                                  # Midpoint shift
  b  <- 1:(N-(B-1))                                   # block indices
  tb <- b + B2 - 1                                    # time in-block
  vb <- t(spec.mtm(ts(1:B), plot = FALSE)$mtm$dpss$v) # DPSS (based on B)
  BS <- list(B = B, B2 = B2, b = b, tb = tb, vb = vb) # store for output
  
  
  ## --------- Choose Signals --------- ##
  if(!is.null(freq)){
    
    # Initialize    
    ind <- 1
    st  <- list()
    
    # st matrix: each column is a signal, rowsums = chord
    for(f in freq){ st[[ind]] = cos(2*pi*Fs[f]*t); ind <- ind + 1 }
    sig <- rowSums(matrix(unlist(st), nrow = N))
    
  } else { st <- NULL; sig <- 0 }
  
  
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
  
  
  # Loop setup
  if(seed > 0){ set.seed(seed) }
  pb <- txtProgressBar(style = 3)
  
  # Initialize arrays
  bc1 <- if(any(p==1)) {array(dim = c(NF,N,M))} else { NULL }
  bc2 <- if(any(p==2)) {array(dim = c(NF,N,M))} else { NULL }
  
  
  
  ## -------------------- SIMULATE: loop begins HERE -------------------- ##
  for(m in 1:M){
    
    xt <- get()
    
    # Sliding window: full window included
    sw <- (B/K) * Mod( sapply(b, function(bb){
      (spec.mtm(ts(xt[(bb):(bb+B-1)]),
                nFFT = (2*NF-1),
                plot = FALSE,
                returnInternals = TRUE)$mtm$eigenCoefs) %*% vb }))^2 
    
    # Endpoints: 1st derivative
    tds.start  <- tds(sw,1)
    tds.end    <- tds(sw,max(b)) 
    
    # Endpoints: 2nd derivative
    stds.start <- stds(sw,1)      
    stds.end   <- stds(sw,max(b))
    
    # Heavy lifting done 
    gc(); setTxtProgressBar(pb, (m-0.5)/M)
    
    
    # --------- 1-BC --------- #
    if(any(p==1)){
      
      # Modify (linear) storing midpoints
      bc1[,tb,m] <- (K/(mval[1]*w0.0)) * sapply(b, function(p){matrix(sw[,p], nrow = NF) %*% mvec[,1]})
      
      # Extrapolate
      bc1[,1:(tb[1]-1),m]   <- sapply((tb[1]-1):1, 
                                      function(h){ bc1[,tb[1],m] - h*tds.start })
      bc1[,(max(tb)+1):N,m] <- sapply(1:(B2-(B2 != B/2)), 
                                      function(h){ bc1[,max(tb),m] + h*tds.end })
    }
    
    
    # --------- 2-BC --------- #
    if(any(p==2)){
      
      # Modify (quadratic) storing midpoints
      bc2[,tb,m] <- (B/d.02) * sapply(b, function(p){ (a(sw,2,p)*w0.2 - a(sw,0,p)*w2.2) })
      
      # Extrapolate
      bc2[,1:(tb[1]-1),m]   <- sapply((tb[1]-1):1, 
                                      function(h){ bc2[,tb[1],m] - h*tds.start - (h^2)*stds.start/2 })
      bc2[,(max(tb)+1):N,m] <- sapply(1:(B2-(B2 != B/2)), 
                                      function(h){ bc2[,max(tb),m] + h*tds.end + (h^2)*stds.end/2 })
    }
    
    # Clean!
    rm(sw); gc()
    
    # Phew!
    setTxtProgressBar(pb, m/M); gc()
  } ## --------------------------- loop ends --------------------------- ##
  
  
  
  # # TRUE values [defunct] ----
  # yf.true  <- if(type != "noise"){AR2.spec(ar2)  } else {NULL}
  # ytf.true <- if(type != "noise"){outer(rep(1,N), yf.true )} else {NULL}
  # 
  # yf2.true  <- if(type == "GUMP"){AR2.spec(ar2.2)} else {NULL}
  # ytf2.true <- if(type == "GUMP"){outer(rep(1,N), yf2.true)} else {NULL}
  # 
  # xtf.true  <- if(type == "AR2"    ){ ytf.true
  # } else if(type == "UMP"    ){ outer(ct1^2,  yf.true)
  # } else if(type == "GUMP"   ){ outer(ct1^2, yf.true) + outer(ct2^2, yf2.true)
  # } else if(type == "AR2_non"){
  #   1e2/(Mod(1 - (a.12)%*%exp(-i*2*pi*Fs) - psi.2%*%exp(-i*2*pi*Fs*2))^2)
  # }
  # xtf1.true <- if(type == "GUMP"){outer(ct1^2,  yf.true)} else {NULL}
  # xtf2.true <- if(type == "GUMP"){outer(ct2^2, yf2.true)} else {NULL}
  # 
  # # Save true vals
  # tru <- list( yf.true, ytf.true, yf2.true, ytf2.true, xtf.true, xtf1.true, xtf2.true)
  # names(tru) <- c("yf",    "ytf",    "yf2",    "ytf2",    "xtf",    "xtf2")
  # ----
  
  # output
  return(list(bc1 = bc1,
              bc2 = bc2,
              # tru = tru, 
              B.vars = BS, 
              sigs = st))
}