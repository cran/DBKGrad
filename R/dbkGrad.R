dbkGrad  <- function(obsqx, omega, ex=NULL, logit=F, bandwidth="FX", h=0.002, s=0.2, cvres="propres", cvh=T, cvs=F, alpha=0.05){
  if (missing(omega)) {omega<-length(obsqx)-1}
  obsqx <- as.vector(obsqx[0:(omega+1)])
  if(!is.null(ex)){ex <- as.vector(ex[0:(omega+1)])}
  if(!any(bandwidth ==c("EX","VC","FX"))) {stop("error: the value specified for the bandwidth parameter is not FX, EX or VC")}  
  if(!any(cvres ==c("propres","res"))) {stop("error: the value specified for the cvres parameter is not propres or res")}
  if(is.null(ex) & any(bandwidth ==c("EX","VC"))) {stop("error: no value specified for the ex parameter")}
  if (bandwidth=="VC"){
    weights <- .vcf(ex,obsqx)
    weights <- weights/sum(weights)
  }
  if (bandwidth=="EX"){
    weights <- sum(ex)/ex
    weights <- weights/max(weights)
  }
  if (bandwidth=="FX"){weights <- rep(1, omega+1)}
  if (logit){
    obsqx <- replace(obsqx, obsqx==0, 1e-6)
    obsqx <- log(obsqx/(1-obsqx))
  }                     
  if (cvh | cvs) {
    require("minpack.lm")
    if (cvh & cvs) {
      nls.out <-  nls.lm(par=c(h,s), cvres=cvres, omega=omega, obsqx=obsqx, weights=weights, fn=.dkCV, lower=c(1e-200,0), upper=c(Inf,1), control = nls.lm.control(maxiter=1000,nprint=2))
      h <- nls.out$par[1]
      s <- nls.out$par[2]
    }
    if (cvh & !cvs) {
      nls.out <-  nls.lm(par=h, s=s, cvres=cvres, omega=omega, obsqx=obsqx, weights=weights, fn=.dkCV, lower=1e-200, upper=Inf, control = nls.lm.control(maxiter=1000,nprint=2)) 
      h <- nls.out$par
    }
    if (!cvh & cvs) {
      nls.out <- nls.lm(par=s, h=h, cvres=cvres, kernel=kernel, omega=omega, obsqx=obsqx, weights=weights, fn=.dkCV, lower=0, upper=1, control = nls.lm.control(maxiter=1000,nprint=2))
      s <- nls.out$par
    }
    cvRSS <- nls.out$rsstrace[nls.out$niter]
    if (is.nan(cvRSS)){
      cat("Warning: for the specified h", ifelse(bandwidth!="FX","and s",""), "cross-validation returned nan; specify different initial values.")}
  }
  else{
    cvRSS <- .dkCV(par=NULL, h=h, s=s, cvres=cvres, omega=omega, obsqx=obsqx, weights=weights)
    if (any(is.nan(cvRSS))){
        cat("Warning: for the specified h", ifelse(bandwidth!="FX","and s,",""),"no smoothing was applied at ages:")
        cat((0:omega)[is.nan(cvRSS)],sep = ",")
    }
    cvRSS <- sum(cvRSS)^2
  }
  K <- .dkern(h=h,s=s, weights=weights, omega=omega)
  kernels  <- K/rowSums(K)
  qxest    <- as.vector(kernels %*% obsqx)
  if (logit) 	{ 
    qxest  <- exp(qxest)/(1+exp(qxest))
    obsqx  <- exp(obsqx)/(1+exp(obsqx))
    }
  result <- list(fitted.values=qxest, h=h, s=s, residuals=qxest-obsqx, prop.residuals=(qxest/obsqx)-1, obsqx=obsqx, omega=omega, ex=ex, kernels=kernels,cvRSS=cvRSS, call=match.call())
  if(!is.null(ex)){
    varqx <- qxest*(1-qxest)/ex
    halfCIlength <- qnorm(1-alpha/2) * (as.vector(kernels^2 %*% varqx))^(1/2)
    upperbound <- sapply(qxest+halfCIlength,function(x)min(x,1))
    lowerbound <- sapply(qxest-halfCIlength,function(x)max(x,1e-100))
    result <- c(result,list(upperbound=upperbound,lowerbound=lowerbound))
  }
  class(result) <- "dbkGrad"
  result
}
as.data.frame.dbkGrad <- function(x, row.names = 0:x$omega, optional = FALSE, ...)
{
  l <- list(obsqx=x$obsqx, fitted.values=x$fitted.values)
  if(!is.null(x$ex)){
    l <- c(l,list(exposed=x$ex, lowerbound=x$lowerbound, upperbound=x$upperbound))
  }
  df <- as.data.frame(l,row.names=row.names, optional=optional, ...)
df
}
print.dbkGrad <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nGraduated Rates:\n")
  print(x$fitted.values)
}

.dkern <- function(h,s, weights, omega) {
  b  <- function(i,j){    
    y <- (x[j]+0.5)/(omega+1)
    a <- (m[i]+0.5)/(bandwidth*(omega+1))+1
    b <- (omega+0.5-m[i])/(bandwidth*(omega+1))+1
    db <- dbeta(x=y,shape1=a,shape2=b,ncp = 0,log = FALSE)
    db
  }
  bandwidth <- rep(h * weights^s, omega+1) #
  x  <- m <- 0:omega
  K  <- array(0,c(omega+1,omega+1),dimnames=list(x,m))
  gridi <- gridj <- 1:(omega+1) 
  K <- outer(gridi,gridj,b)
  K
}
.dkCV <- function(par, h, s, bandwidth, omega, obsqx, cvres, weights) {
  if (missing(h) & missing(s)){
    h <- par[1]
    s <- par[2]
  }
  if (!missing(h) & missing(s)){
        s <- par[1]
      }    
  if (missing(h) & !missing(s)){
        h <- par[1]
  }    
  K <- .dkern(h=h, s=s, weights=weights, omega=omega)
  Kd <- K - diag(diag(K))
  Kdrop <- Kd/rowSums(Kd)
  obsqxremove <- Kdrop %*% obsqx 
  if (cvres=="propres") { CV = obsqxremove/obsqx-1 }
  else {CV = obsqxremove-obsqx}
  CV
}

.vcf <- function(n,p){
  sqrt((1-p)/(n*p))
}
plot.dbkGrad <- function(x, plottype="obsfit", CI=T,...)
{
  if(!any(plottype==c("obsfit","fitted", "observed","histres","histpropres", "exposed" ))){
  stop("error: the value specified for the bandwidth parameter is not obsfit, fitted, observed, histres, histpropres or exposed")}
  eval(parse(text=paste0(".",plottype,".plot(x=x,...)")))
  if(CI & any(plottype==c("obsfit","fitted")& !is.null(x$ex))){.ci.plot(x=x,...)}
}

.histres.plot<- function(x, ...)
{
 # plot(density(residuals(x)),main="",xlab="residuals",...)
  hist(residuals(x),main="",xlab="residuals",...)
}
.histpropres.plot<- function(x, ...)
{
  # plot(density(residuals(x)),main="",xlab="residuals",...)
  hist(x$prop.residuals,main="",xlab="residuals",...)
}
.ci.plot <-function(x,...)
{
    segments(x0=0:x$omega,y0=log(x$fitted.values),x1=0:x$omega,y1=log(x$upperbound),col="black")
    segments(x0=0:x$omega,y0=log(x$fitted.values),x1=0:x$omega,y1=log(x$lowerbound),col="black")
}
.obsfit.plot<- function(x,...)
{
  qxmat <- cbind(log(x$fitted.values),log(x$obsqx))
  grid  <- seq(min(qxmat),max(qxmat),length=5)
  par(mai=c(0.9,0.9,0.1,0.02))
  par(cex.lab = 1.4)  
  par(cex.axis = 0.9) 
  par(las = 3)
  matplot(x=0:x$omega,y=qxmat,type="p",pch=c(16, 1),xlab="Age",ylab="mortality rates (log scale)",col=(c("black","blue")),axes = FALSE, ...)
  axis(1,at = seq(0,x$omega,by=5),labels = round(seq(0,x$omega,by=5)))
  axis(2,at = grid, labels = round(exp(grid),digits=5))
  box(col = "black")
  legend("bottomright",legend = c("graduated","observed"), pch=c(16, 1), col=(c("black","blue")),inset=0.02,text.width=strwidth("grad"))		
}
.observed.plot<- function(x,...)
{
  grid <- seq(min(log(x$obsqx)),max(log(x$obsqx)),length=5)
  par(mai=c(0.9,0.9,0.1,0.02))
  par(cex.lab = 1.4)  
  par(cex.axis = 0.9) 
  par(las = 3)
  plot(x=0:x$omega,y=log(x$obsqx),type="p",pch=1,xlab="Age",ylab="observed mortality rates (log scale)",col="blue",axes = FALSE, ...)
  axis(1,at = seq(0,x$omega,by=5),labels = round(seq(0,x$omega,by=5)))
  axis(2,at = grid, labels = round(exp(grid),digits=5))
  box(col = "black")		
}
.fitted.plot<- function(x,...)
{
  grid <- seq(min(log(x$fitted.values)),max(log(x$fitted.values)),length=5)
  par(mai=c(0.9,0.9,0.1,0.02))
  par(cex.lab = 1.4)  
  par(cex.axis = 0.9) 
  par(las = 3)
  plot(x=0:x$omega,y=log(x$fitted.values),type="p",pch=16,xlab="Age",ylab="graduated mortality rates (log scale)",col="black",axes = FALSE, ...)
  axis(1,at = seq(0,x$omega,by=5),labels = round(seq(0,x$omega,by=5)))
  axis(2,at = grid, labels = round(exp(grid),digits=5))
  box(col = "black")		
}
.exposed.plot<- function(x,...)
{
  grid <- seq(0,max(x$ex),length=5)
  par(mai=c(0.9,0.9,0.1,0.02))
  par(cex.lab = 1.4)  
  par(cex.axis = 0.9) 
  par(las = 3)
  plot(x=0:x$omega, y=x$ex, xlab="Age", ylab="exposed to risk", type="h", lwd=2, axes = FALSE, ...)
  axis(1,at = seq(0,x$omega,by=5),labels = round(seq(0,x$omega,by=5)))
  axis(2,at = grid, labels = round(grid,digits=5))
  box(col = "black")		
}