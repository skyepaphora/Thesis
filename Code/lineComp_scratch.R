






# ---------- The random version? ----
gest.rand <- function(sgram, n){
  N  <- length(sgram[,1])
  NF <- length(sgram[1,])
  l <- sample(1:NF, n)
  A <- outer(1:N, 1:N, function(j,k){
    rowSums(sgram[,l])[j]/rowSums(sgram[,l])[k]
  })
  G <- Mod(eigen(A/N)$vectors[,1])
  return(G); gc()
}

set.seed(11)
n <- 50
L <- 50
out <- matrix(0, N, L)

plot.true()
for(m in 1:L){
  out[,m] <- gest.rand(sgram,n)
  lines(t, normy(out[,m]), col = "#0000FF33")  
  setTxtProgressBar(pb, m/L)
}
gc()


lines(t, old.g, col = "goldenrod1", lwd = 2)
lines(t, svd.g, col = "red", lwd = 2)
lines(t, normy(rowMeans(out)), col = "cyan", lwd = 2, lty = 2)
legend("topright",
       c("true","ours (new)","E[ours (new)]", "ours (old)", "SVD"),
       col = c("black","blue","cyan","goldenrod","red"),
       lty = c(1,1,2,1,1), cex = 0.8, lwd = c(1,1,2,2,2))

mtext(paste("random column selection  | ",
            L, "sets of", n, "cols"))

MSE.rand <- mean(abs(normy(rowMeans(out)) - normy(gt.1)))
# MSE.n100 <- mean(abs(normy(rowMeans(out)) - normy(gt.1)))
legend("bottomleft", paste("n =",n,"; MSE =", MSE.rand))




# 0. tester
## sgram <- t(pbc.test$bc1[,,1]) # t(pbc.go)

# 1. Lifted Cosine ----
set.seed(11)
n <- 200
O <- 200

omega <- f[O]
c.t   <- cos(2*pi*omega*t) + 1.5 + ct$grow$line
phii  <- c(0.75,-0.5)

pbc.mod <- pbc(ct.1 = c.t, s.t = s.t, sigs.mod = 1, M = 1, p = 1, seed = 1)
Sx      <- t(pbc.mod$bc1[,,1])
xt      <- pbc.mod$xt.last$xt

g.sunhat <- jester(Sx, c.t, n, plot = FALSE)
g.hat <- g.sunhat$mean
c.hat <- sqrt(g.hat - (min(g.hat)-1))
y.hat <- xt/c.hat

x.mtm <- spec.mtm(xt   , Ftest = TRUE, plot = FALSE)
y.mtm <- spec.mtm(y.hat, Ftest = TRUE, plot = FALSE)

splot(f, x.mtm$mtm$Ftest,  type = 'h', lwd = 3, skor = FALSE)
lines(f, y.mtm$mtm$Ftest, col = "green", lty = 3, lwd = 2, type = "h")
mtext(paste("n =", n, " |  Fourier freq =", O), line = 2)
gc()

# 2. smoo scales ----
gg <- gest(Sx)
ss <- yest(Sx)

range(outer(c.t^2, yf))
range(outer(gg, ss))

summary(gg) *N/12
summary(c.t^2)
summary(ss) *NF/12
summary(yf)

# 3. Apply Glen's code to UMP + arbitrary UMP ----

set.seed(11)
phii  <- c(0.75,-0.5)

ct.1 <- ct$grow$sqrt
ct.2 <- ct$decay$line
c1st <- normy(ct.1)
c2st <- normy(ct.2)

pbc.mod <- pbc(ct.1 = ct.1, ct.2 = ct.2,
               s.t = s.t, sigs.mod = 1,
               M = 1, p = 1, seed = 1)
gump2 <- t(pbc.mod$bc1[,,1])

j <- sample(1:NF,2)
S1 <- JADE(cbind(gump2[,min(j)],gump2[,max(j)]))$S

par(mfrow = c(2,1))
plot(S1[,1],type="l",col="blue", main = paste("sample = ", min(j)))
lines(c1st,col="red")
plot(S1[,2],type="l",col="blue", main = paste("sample = ", max(j)))
lines(c2st,col="red")




