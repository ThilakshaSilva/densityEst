#' @title Plotting density estimation using logspline approach
#' @description  It nonparametrically estimates a time series of density functions using logspline approach considering the time ordering of the densities.
#' @param final A list of outputs returns from densityEst function
#' @export
plots <- function(final)
{
  y <- final$y
  den <- final$z
  optp <- final$optp
  knotsy <- final$optimalknots
  ymargin <- final$ymargin
  xmargin <- final$xmargin
  nymargin <- final$nymargin
  nxmargin <- final$nxmargin


  logden <- log(den)
  logknotsy <- log(knotsy)
  logymargin <- log(ymargin)
  den1 <- sapply(1:nxmargin, function(X) den[,X]*ymargin)
  logden1 <- log(den1)
  ly <- sapply(1:nxmargin, function(X) length(y[[X]]))


  graphics::matplot(ymargin[1:nxmargin],den[1:nxmargin,],type="l",col=grDevices::rainbow(nxmargin,start=0,end=4/6),lty=1,ylab="Density",main=paste("Density functions (",xmargin[1],"-",xmargin[nxmargin],")"))
  graphics::legend("topright",c(paste(xmargin[1]),paste(xmargin[nxmargin])),col=c("red","blue"),lty=1)
  graphics::rug(knotsy,lwd=2,col="red")


  a <- list()
  a$x <- xmargin
  a$y <- ymargin[1:nxmargin]
  a$z <- t(den[1:nxmargin,])
  graphics::par(mfrow=c(1,1))
  hdrcde::plot.cde(a,expand=0.3,xlab="Year",col="dodgerblue")
}


