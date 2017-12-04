"basis" <- function(y,t,d,e=TRUE)
{
  terms <- matrix(NA,length(y),length(t))
  if (e==FALSE)
  { cubic <- function(y,t) { (t>y)*(t-y)*(t-y)*(t-y) } }
  else { cubic <- function(y,t) { (y>t)*(y-t)*(y-t)*(y-t) } }
  terms[,1:length(t)] <- outer(y,t,FUN="cubic")
  bas <- terms%*%d
  max <- max(bas)   #I divide the first and last two basis function by their maximum values to scale them to have values between [0,1]
  bas <- bas/max
  return(list(bas=bas,max=max))
}


### Computing the maximum values of first and last two bases functions
"maxbasiscal" <- function(y,knots)
{
  K <- length(knots)
  max <- numeric(length=3)


  dterms1 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[2]-t[1])
    return(c(-d2-1,d2,1))
  }

  t <- c(knots[1],knots[2],knots[3])
  max[1] <- basis(y,t,d=dterms1(t),e=FALSE)$max

  t <- c(t,knots[4])
  q1 <- t[4]-t[1]
  q2 <- t[2]-t[1]
  d <- numeric(length=4)
  d[2] <- q1*(t[4]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[3]-t[1])
  d[4] <- 1
  d[1] <- -d[2]-d[3]-d[4]
  max[4] <- basis(y,t,d,e=FALSE)$max


  dterms2 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[3]-t[2])
    return(c(1,d2,-d2-1))
  }

  t <- c(knots[K-2],knots[K-1],knots[K])
  max[2] <- basis(y,t,d=dterms2(t))$max

  t <- c(knots[K-3],t)
  q1 <- t[4]-t[1]
  q2 <- t[4]-t[2]
  d[1] <- 1
  d[2] <- q1*(t[1]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[4]-t[3])
  d[4] <- -d[1]-d[2]-d[3]
  max[3] <- basis(y,t,d)$max

  return(max)
}


### First derivative of bases functions
"firstderbasis" <- function(y,t,d,e=TRUE)
{
  terms <- matrix(NA,length(y),length(t))
  if (e==FALSE)
  { cubic <- function(y,t) { (t>y)*(t-y)*(t-y) } }
  else { cubic <- function(y,t) { (y>t)*(y-t)*(y-t) } }
  terms[,1:length(t)] <- outer(y,t,FUN="cubic")
  bas <- terms%*%d
  return(bas)
}


"firstder" <- function(y,knots)
{
  y <- unlist(y)
  y <- sort(y)
  K <- length(knots)
  max <- maxbasiscal(y,knots)


  dbas <- matrix(NA,length(y),K-1)
  ### Following basis functions are exactly same - gives 12 basis functions for 11 knots
  #dbas[,2:(K-3)] <- dbasisgen(x,degree=3,interior.knots=knots[2:(K-1)],Boundary.knots=c(knots[1],knots[K]))[,3:(K+1-3)]
  dbas[,3:(K-3)] <- crs::gsl.bs(y,intercept=FALSE,knots=knots,nbreak=K,deriv=1)[,4:(K+1-3)]


  dterms1 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[2]-t[1])
    return(c(-3*(-d2-1),-3*d2,-3))
  }

  t <- c(knots[1],knots[2],knots[3])
  dbas[,1] <- firstderbasis(y,t,d=dterms1(t),e=FALSE)/max[1]   #Generate B_1 B-spline function

  t <- c(t,knots[4])
  q1 <- t[4]-t[1]
  q2 <- t[2]-t[1]
  d <- numeric(length=4)
  d[2] <- q1*(t[4]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[3]-t[1])
  d[4] <- 1
  d[1] <- -d[2]-d[3]-d[4]
  dbas[,2] <- firstderbasis(y,t,(-3)*d,e=FALSE)/max[4]         #Generate B_2 B-spline function


  dterms2 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[3]-t[2])
    return(c(3,3*d2,3*(-d2-1)))
  }

  t <- c(knots[K-2],knots[K-1],knots[K])
  dbas[,K-1] <- firstderbasis(y,t,d=dterms2(t))/max[2]         #Generate B_{K-1} B-spline function

  t <- c(knots[K-3],knots[K-2],knots[K-1],knots[K])
  q1 <- t[4]-t[1]
  q2 <- t[4]-t[2]
  d[1] <- 1
  d[2] <- q1*(t[1]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[4]-t[3])
  d[4] <- (-d[1]-d[2]-d[3])
  dbas[,K-2] <- firstderbasis(y,t,3*d)/max[3]                    #Generate B_{K-2} B-spline function

  return(dbas)
}


### Second derivative of bases functions
"secondderbasis" <- function(y,t,d,e=TRUE)
{
  terms <- matrix(NA,length(y),length(t))
  if (e==FALSE)
  { cubic <- function(y,t) { (t>y)*(t-y) } }
  else { cubic <- function(y,t) { (y>t)*(y-t) } }
  terms[,1:length(t)] <- outer(y,t,FUN="cubic")
  bas <- terms%*%d
  return(bas)
}


"secondder" <- function(y,knots)
{
  y <- unlist(y)
  y <- sort(y)
  K <- length(knots)
  max <- maxbasiscal(y,knots)


  d2bas <- matrix(NA,length(y),K-1)
  ### Following basis functions are exactly same - gives 12 basis functions for 11 knots
  #d2bas[,2:(K-3)] <- d2basisgen(y,degree=3,interior.knots=knots[2:(K-1)],Boundary.knots=c(knots[1],knots[K]))[,3:(K+1-3)]
  d2bas[,3:(K-3)] <- crs::gsl.bs(y,intercept=FALSE,knots=knots,nbreak=K,deriv=2)[,4:(K+1-3)]


  dterms1 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[2]-t[1])
    return(c(6*(-d2-1),6*d2,6))
  }

  t <- c(knots[1],knots[2],knots[3])
  d2bas[,1] <- secondderbasis(y,t,d=dterms1(t),e=FALSE)/max[1]         #Generate B_1 B-spline function

  t <- c(t,knots[4])
  q1 <- t[4]-t[1]
  q2 <- t[2]-t[1]
  d <- numeric(length=4)
  d[2] <- q1*(t[4]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[3]-t[1])
  d[4] <- 1
  d[1] <- -d[2]-d[3]-d[4]
  d2bas[,2] <- secondderbasis(y,t,6*d,e=FALSE)/max[4]                  #Generate B_2 B-spline function


  dterms2 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[3]-t[2])
    return(c(6,6*d2,6*(-d2-1)))
  }

  t <- c(knots[K-2],knots[K-1],knots[K])
  d2bas[,K-1] <- secondderbasis(y,t,d=dterms2(t))/max[2]               #Generate B_{K-1} B-spline function

  t <- c(knots[K-3],t)
  q1 <- t[4]-t[1]
  q2 <- t[4]-t[2]
  d[1] <- 1
  d[2] <- q1*(t[1]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[4]-t[3])
  d[4] <- -1-d[2]-d[3]
  d2bas[,K-2] <- secondderbasis(y,t,6*d)/max[3]                        #Generate B_{K-2} B-spline function

  return(d2bas)
}


### Basis calculation
"basiscaly"  <- function(y,knotsy)
{
  y <- unlist(y)
  y <- sort(y)
  K <- length(knotsy)

  z <- matrix(NA,length(y),K-1)


  ### Following basis functions are exactly same - gives 7 basis functions for 11 knots
  ### Remove first two bases functions and last three bases functions
  if((K-2) < 3)       #K < 5
    stop("At least 5 knots are needed")
  if((K-2) == 3)      #K = 5
  {
    #z[,2] <- basisgen(y,intercept=FALSE,interior.knots=knots[2:(K-1)],Boundary.knots=c(knots[1],knots[K]))[,3]
    z[,2] <- crs::gsl.bs(y,intercept=FALSE,knots=knotsy,nbreak=K)[,3]
  }
  if((K-2) != 3)      #K != 5
  {
    #z[,2:(K-3)] <- basisgen(y,intercept=FALSE,interior.knots=knots[2:(K-1)],Boundary.knots=c(knots[1],knots[K]))[,3:(K-2)]
    z[,3:(K-3)] <- crs::gsl.bs(y,intercept=FALSE,knots=knotsy,nbreak=K)[,4:(K-2)]
  }


  dterms1 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[2]-t[1])
    return(c(-d2-1,d2,1))
  }

  t <- c(knotsy[1],knotsy[2],knotsy[3])
  z[,1] <-  basis(y,t,d=dterms1(t),e=FALSE)$bas    #Generate B_1 B-spline function

  t <- c(t,knotsy[4])
  q1 <- t[4]-t[1]
  q2 <- t[2]-t[1]
  d <- numeric(length=4)
  d[2] <- q1*(t[4]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[3]-t[1])
  d[4] <- 1
  d[1] <- -d[2]-d[3]-d[4]
  z[,2] <- basis(y,t,d,e=FALSE)$bas                #Generate B_2 B-spline function


  dterms2 <- function(t)
  {
    d2 <- (t[1]-t[3])/(t[3]-t[2])
    return(c(1,d2,-d2-1))
  }

  t <- c(knotsy[K-2],knotsy[K-1],knotsy[K])
  z[,K-1] <- basis(y,t,d=dterms2(t))$bas           #Generate B_{K-1} B-spline function

  t <- c(knotsy[K-3],knotsy[K-2],knotsy[K-1],knotsy[K])
  q1 <- t[4]-t[1]
  q2 <- t[4]-t[2]
  d[1] <- 1
  d[2] <- q1*(t[1]-t[3])/(q2*(t[3]-t[2]))
  d[3] <- (-q1-q2*d[2])/(t[4]-t[3])
  d[4] <- -1-d[2]-d[3]
  z[,K-2] <- basis(y,t,d)$bas                      #Generate B_{K-2} B-spline function

  return(z)
}


### Maximum likelihood - maximize log density values
"lik" <- function(intp,y,basisy,ymargin,delta)
{
  zb <- basisy%*%intp
  cb <- log(sum(exp(zb))*delta)
  z <- stats::approx(ymargin,zb-cb,y)$y        #log density of log income
  return(sum(z))
}


### Choosing optimal parameters
"choosepar" <- function(b1,intp,y,basisy,ymargin,delta,K)
{

  b <- optimx::optimx(intp,lik,method="nlminb",control=list(maximize=TRUE),upper=c(b1,rep(Inf,K-3),-0.1),y=y,basisy=basisy,ymargin=ymargin,delta=delta) #Calls lik function and after repeating provide maximum likelihood
  return(list(par=as.vector(stats::coef(b)),value=b$value))     #Maximixe the optimx function by fnscale=-1 and returns 'length(knots)+1' number of best set of coefficients
}


### Obtaining mid knots
"allknots" <- function(y,K,lbound,ubound)
{
  #Modified Chebyshev knots
  #Knots are calculated in original scale
  data <- exp(sort(as.vector(unlist(y))))
  point1 <- 1
  point2 <- length(data)
  intk <- sapply(1:K, function(x) (point1+point2)/2 + ((point2-point1)/2)*cos((x-1)*pi/(K-1)))
  intk <- data[rev(intk)]

  intk[1] <- exp(lbound)
  intk[K] <- exp(ubound)

  midk <- numeric(length=(K-3))
  midk[1:(K-3)] <- (intk[2:(K-2)] + intk[3:(K-1)])/2

  ##########################################
  #fmidk <- (intk[1] + intk[2])/2

  #Knots are calculated again in original scale
  ###The spike of the second basis is always between first and second inner knots.
  ###Find the first mid knot which is between first and second inner knots. In the graph, this mid knot is in usually after the spike of the second basis.
  ###I find new 12 inner knots, after the first mid knot.
  ###First inner knot is between old first mid knot and old first inner knot.
  ###So new first inner knot is in usually after the spike of t                                                               he second basis. Now it is assumed that no unusual behavior arise from this new first inner knot as it is after the spike of the second basis.
  #point <- (1:length(data))[data>fmidk][1]
  #int <- numeric(length=nknots)
  #for(i in 1:nknots)
  #{
  #  int[i] <- (point+length(data))/2 + ((length(data)-point)/2)*cos((2*i-1)*pi/(2*nknots))
  #}
  #intk <- data[rev(int)]
  #midk <- numeric(length=nknots)
  #midk[1:(nknots-1)] <- (intk[1:(nknots-1)] + intk[2:nknots])/2
  #midk[nknots] <- (intk[nknots]+data[length(data)])/2

  return(list(intk=log(intk),midk=log(midk)))
}


### Starting values - Initial parameters
"intpar" <- function(y,knotsy)
{
  #matplot(sort(unlist(y)),basiscaly(y,knotsy),type="l")
  #rug(knotsy)
  fdif <- firstder(y,knotsy)               #n*(length(knotsy)-1)
  #matplot(sort(unlist(y)),fdif,type="l")
  #rug(knotsy)
  sdif <- secondder(y,knotsy)              #n*(length(knotsy)-1)
  #matplot(sort(unlist(y)),sdif,type="l")
  #rug(knotsy)

  left <- t(fdif)%*%fdif
  right <- -colSums(sdif)
  intp <- solve(left,right)

  b1 <- 1/fdif[1,1]
  if (intp[1] > b1)
    intp[1] <- b1
  K <- length(knotsy)
  if (intp[K-1] > -0.5)
    intp[K-1] <- -0.5

  return(list(intp=intp,b1=b1))
}


###BIC for each dataset with different initial parameters
#If you want to export your complete global environment then just write .export=ls(envir=globalenv()) and you will have it for better or worse.
"biccal" <- function(b1,intp,y,knotsy,ymargin,delta,nxmargin,nmidk,n=2)
{
  doParallel::registerDoParallel(n)

  K <- length(knotsy[1,])
  i = j = 0

  `%:%` <- foreach::`%:%`
  `%dopar%` <- foreach::`%dopar%`
  bic1 <- foreach::foreach(i=1:nmidk) %:%
    foreach::foreach(j=1:nxmargin , .export=c("basis", "basiscaly", "lik","choosepar")) %dopar%
    choosepar(b1[i],intp[i,],y[[j]],basiscaly(ymargin,knotsy[i,]),ymargin,delta,K)
  bic <- sapply(1:nmidk, function(X1) sapply(1:nxmargin, function(X2) -2*bic1[[X1]][[X2]]$value + (K-1)*log(length(y[[X2]]))))
  optp <- lapply(1:nmidk, function(X1) t(sapply(1:nxmargin, function(X2) bic1[[X1]][[X2]]$par)))
  #   bic <- matrix(NA,nxmargin,nmidk)
  #   optp <- vector(length=nmidk,mode="list")
  #   for(j in 1:nmidk)
  #   {
  #     basisy <- basiscaly(ymargin,knotsy[j,])
  #     bic1 <- foreach(i = 1:nxmargin, .export=c("choosepar","lik")) %dopar%
  #       choosepar(b1[j],intp[j,],y[[i]],knotsy[j,],basisy,ymargin,delta,K)
  #
  #     bic[,j] <- sapply(1:nxmargin, function(X) -2*bic1[[X]]$value + K*log(length(y[[X]])))
  #     optp[[j]] <- t(sapply(1:nxmargin, function(X) bic1[[X]]$par))
  #   }

  return(list(bic=bic,knotsy=knotsy,optp=optp))
}


### Knot addition
"addition" <- function(y,K,lbound,ubound,ymargin,nxmargin,delta,n=2)
{
  doParallel::registerDoParallel(n)

  ### Initial knots
  knots <- allknots(y,K,lbound,ubound)
  initialknots <- knots$intk
  knotsy <- initialknots


  ### Starting values
  i = 0
  intpp <- intpar(y,knotsy)
  intp <- intpp$intp
  b1 <- intpp$b1                           #b1 <- 1/firstder(y,knotsy)[1,1]
  `%dopar%` <- foreach::`%dopar%`
  bic1 <- foreach::foreach(i = 1:nxmargin, .export=c("basis", "basiscaly", "lik", "choosepar")) %dopar%
    choosepar(b1,intp,y[[i]],basiscaly(ymargin,knotsy),ymargin,delta,K)
  bic <- sapply(1:nxmargin, function(X) -2*bic1[[X]]$value + (K-1)*log(length(y[[X]])))
  optp <- t(sapply(1:nxmargin, function(X) bic1[[X]]$par))
  oldbicfinal <- mean(bic)


  ### Addition of the first mid knot
  midk <- knots$midk
  nmidk <- length(midk)
  newknots <- t(sapply(1:nmidk, function(x) sort(c(knotsy,midk[x]))))

  intpp <- lapply(1:nmidk, function(X) intpar(y,newknots[X,]))
  intp <- t(sapply(1:nmidk, function(X) intpp[[X]]$intp))
  b1 <- sapply(1:nmidk, function(X) intpp[[X]]$b1)
  #   intp <- vector(length=nmidk,mode="list")
  #   b1 <- numeric(length=nmidk)
  #   for(i in 1:nmidk)
  #   {
  #     intpp <- intpar(y,newknots[i,])
  #     intp[[i]] <- intpp$intp
  #     b1[i] <- intpp$b1
  #   }
  bic <- biccal(b1,intp,y,newknots,ymargin,delta,nxmargin,nmidk)
  bicfinal <- colMeans(bic$bic)


  ### Addition of other knots
  kp <- c(1:nmidk)                                #1 2 3 4 5 6 7 8
  avgs <- round(length(unlist(y))/length(y))
  Kmax <- round(min(3*avgs^(1/5),avgs/4,30))      #Maximum number of knots = 18
  dif <- Kmax-K                                   #18-11=7
  dif <- min(dif,max(kp))                         #7
  #When e=1, number of knots is 12 (i.e. the 1st mid knot is added)
  #When e=7, number of knots is 18 (i.e. the 7th mid knot is added)
  for (e in 1:dif)
  {
    if(min(bicfinal) <= oldbicfinal)
    {
      oldbicfinal <- min(bicfinal)
      point <- match(oldbicfinal,bicfinal)
      kp <- kp[-point]
      midk <- midk[-point]
      nmidk <- length(midk)

      knotsy <- bic$knots[point,]
      newknots <- t(sapply(1:nmidk,function(X) sort(c(knotsy,midk[X]))))

      optp <- bic$optp[[point]]
      intpp <- lapply(1:nmidk, function(X) intpar(y,newknots[X,]))
      intp <- t(sapply(1:nmidk, function(X) intpp[[X]]$intp))
      b1 <- sapply(1:nmidk, function(X) intpp[[X]]$b1)

      bic <- biccal(b1,intp,y,newknots,ymargin,delta,nxmargin,nmidk)
      bicfinal <- colMeans(bic$bic)
    }
    else{break}
  }

  return(list(knotsy=knotsy,initialknots=initialknots,optp=optp,oldbicfinal=oldbicfinal))
}


### Knot deletion
"deletion" <- function(y,K,lbound,ubound,ymargin,nxmargin,delta)
{
  add <- addition(y,K,lbound,ubound,ymargin,nxmargin,delta)
  initialknots <- add$initialknots
  additionknots <- add$knotsy
  knotsy <- add$knotsy
  optp <- add$optp
  oldbicfinal <- add$oldbicfinal


  K <- length(knotsy)
  nmidk <- K-2
  newknots <- t(sapply(1:nmidk,function(x) knotsy[-(1+x)]))
  intpp <- lapply(1:nmidk, function(X) intpar(y,newknots[X,]))
  intp <- t(sapply(1:nmidk, function(X) intpp[[X]]$intp))
  b1 <- sapply(1:nmidk, function(X) intpp[[X]]$b1)
  bic <- biccal(b1,intp,y,newknots,ymargin,delta,nxmargin,nmidk)
  bicfinal <- colMeans(bic$bic)


  #At least 4 knots are needed in generating bases functions. So if there are at least 5 knots left, remove a knot.
  deldif <- K-6
  for(e in 1:deldif)
  {
    if(e!=deldif)
    {
      if(min(bicfinal) <= oldbicfinal)
      {
        oldbicfinal <- min(bicfinal)
        point <- match(oldbicfinal,bicfinal)
        optp <- bic$optp[[point]]

        nmidk <- nmidk-1

        knotsy <- bic$knotsy[point,]
        newknots <- t(sapply(1:nmidk,function(X) knotsy[-(1+X)]))

        intpp <- lapply(1:nmidk, function(X) intpar(y,newknots[X,]))
        intp <- t(sapply(1:nmidk, function(X) intpp[[X]]$intp))
        b1 <- sapply(1:nmidk, function(X) intpp[[X]]$b1)
        bic <- biccal(b1,intp,y,newknots,ymargin,delta,nxmargin,nmidk)
        bicfinal <- colMeans(bic$bic)
      }
      else {break}
    }else if(e==deldif)
    {
      if(min(bicfinal) <= oldbicfinal)
      {
        oldbicfinal <- min(bicfinal)
        point <- match(oldbicfinal,bicfinal)
        optp <- bic$optp[[point]]
        knotsy <- bic$knotsy[point,]
      }
      else {break}
    }
  }

  return(list(optp=optp,knotsy=knotsy,optimalknots=exp(knotsy),additionknots=exp(additionknots),initialknots=exp(initialknots)))
}


#' @title Density estimation using logspline approach
#' @description  It nonparametrically estimates a time series of density functions using logspline approach considering the time ordering of the densities.
#' @param y a time series of data
#' @param x a sequence of numbers from 1 to n where n is the number of time series
#' @param ymargin data interval in the y-direction
#' @param lbound lower bound for the support for the density
#' @param ubound upper bound for the support for the density
#' @param K number of initial number of knots
#' @param nymargin desired length of the seqence
#' @param log Does the logarithmic transformation apply to data? Default is TRUE.
#' @param graphics A logical argument for plots
#' @export
#' @examples
#' l <- 20
#' n <- 2000
#' z <- seq(1000,5000,length.out=n)
#' mean <- seq(mean(z)-400,mean(z)+200,length.out=l)
#' sd <- seq(200,400,length.out=l)
#' y <- lapply(1:l, function(X) rnorm(z,mean[X],sd[X]))
#' x <- c(1:l)
#' densityEst(y,x)

"densityEst" <- function(y,x,ymargin,lbound,ubound,K,nymargin=1000,log=TRUE,graphics=TRUE)
{
  if(missing(y))
    stop("You must provide y")

  min <- min(unlist(y))
  max <- max(unlist(y))
  if (log==TRUE)
  {
    if(missing(ymargin))
    {
      #lbound is the lower bound for the support for the density.
      if(!missing(lbound)){
        if (lbound <= 0)
          stop("lbound can only be a positive value")
        lbound <- log(lbound)
      }else{
        range <- (max - min)/4
        if ((min - range) <= 0)
          lbound <- log(min)
        if ((min - range) > 0)
          lbound <- log(min - range)
      }

      #ubound is the upper bound for the support for the density.
      if(!missing(ubound)){
        if (ubound <= 0)
          stop("ubound can olny be a positive value")
        if (ubound < max(unlist(y)))
          stop("ubound cannot be less than the maximum of all the data")
        ubound <- log(ubound)
      }else
        ubound <- log(max + range)

      ymargin <- seq(lbound,ubound,length.out=nymargin)

    }else{
      #lbound is the lower bound for the support for the density.
      if(!missing(lbound)){
        if (lbound <= 0)
          stop("lbound can only be a positive value")
        lbound <- log(lbound)
      }else{
        range <- (max - min)/4
        if ((min - range) <= 0)
          lbound <- log(min)
        if ((min - range) > 0)
          lbound <- log(min - range)
      }

      #ubound is the upper bound for the support for the density.
      if(!missing(ubound)){
        if (ubound <= 0)
          stop("ubound can olny be a positive value")
        if (ubound < max(unlist(y)))
          stop("ubound cannot be less than the maximum of all the data")
        ubound <- log(ubound)
      }else
        ubound <- log(max + range)
    }

    if(!missing(x)){
      nxmargin <- length(x)
    }

    #Log y values
    if(!missing(y)){
      y <- lapply(1:nxmargin, function(X) log(sort(y[[X]])))
    }
  }else{
    #lbound is the lower bound for the support for the density.
    if(missing(lbound)){
      range <- (max - min)/4
      lbound <- min - range
    }

    #ubound is the upper bound for the support for the density.
    if(missing(ubound)){
      range <- (max - min)/4
      ubound <- max + range
    }

    ymargin <- seq(lbound,ubound,length.out=nymargin)
  }


  #K is the number of initial number of knots.
  if(missing(K)){
    avgssize <- round(length(unlist(y))/length(y))
    K <- round(min(2*avgssize^(1/5),avgssize/4,25))     #K=14(12 inner initial knots + 2 bounday knots)
  }


  #   range <- (max(unlist(y)) - min(unlist(y)))/4
  #   if ((min(unlist(y)) - range) <= 0){
  #     lbound <- log(min(unlist(y)))
  #   }
  #   if ((min(unlist(y)) - range) > 0){
  #     lbound <- log(min(unlist(y)) - range)
  #   }
  #   ubound <- log(max(unlist(y)) + range)
  #   nxmargin <- length(y)
  #   y <- lapply(1:nxmargin, function(X) log(y[[X]]))
  #   avgssize <- round(length(unlist(y))/length(y))
  #   K <- round(min(2*avgssize^(1/5),avgssize/4,25))
  #   nymargin <- 1000
  #   ymargin <- seq(lbound,ubound,length.out=nymargin)
  delta <- ymargin[2]-ymargin[1]
  xmargin <- x
  #####################################################


  #Knot deletion
  del <- deletion(y,K,lbound,ubound,ymargin,nxmargin,delta)
  optp <- del$optp
  knotsy <- del$knotsy
  optimalknots <- del$optimalknots
  additionknots <- del$additionknots
  initialknots <- del$initialknots


  #Density estimation
  basisy <- basiscaly(ymargin,knotsy)   #dim(basisy)=nymargin*(length(knotsy)-1), dim(optp)=nxmargin*(length(knotsy)-1)
  expymargin <- exp(ymargin)

  zb <- sapply(1:nxmargin, function(X) basisy%*%optp[X,])             #log density of log income, dim(zb)=nymargin*nxmargin
  cb <- sapply(1:nxmargin, function(X) log(sum(exp(zb[,X]))*delta))
  z <- sapply(1:nxmargin, function(X) exp(zb[,X]-cb[X])/expymargin)   #density of log income / expymargin = density of income
  int <- sapply(1:nxmargin, function(X) sum(((z[1:(nymargin-1),X]+z[2:nymargin,X])/2)*(expymargin[2:nymargin]-expymargin[1:(nymargin-1)])))
  z <- sapply(1:nxmargin, function(X) z[,X]/int[X])
  #   z <- matrix(NA,nxmargin,nymargin)
  #   for (i in 1:nxmargin)
  #   {
  #     zb <- basisy%*%optp[i,]                                          #log density of log income
  #     cb <- log(sum(exp(zb))*delta)
  #     z[i,] <- exp(zb-cb)                                              #density of log income
  #     z[i,] <- z[i,]/expymargin                                        #density of income
  #     z[i,] <- z[i,]/sum(((z[i,1:(nymargin-1)]+z[i,2:nymargin])/2)*(expymargin[2:nymargin]-expymargin[1:(nymargin-1)]))
  #   }


  #Original y values
  y <- lapply(1:nxmargin, function(X) exp(y[[X]]))

  density <- list()
  density$y <- y
  density$z <- z
  density$optp <- optp
  density$basisy <- basisy
  density$optimalknots <- optimalknots
  density$additionknots <- additionknots
  density$initialknots <- initialknots
  density$ymargin <- expymargin
  density$xmargin <- xmargin
  density$nymargin <- nymargin
  density$nxmargin <- nxmargin
  density$lowerbound <- exp(lbound)
  density$upperbound <- exp(ubound)
  density$delta <- delta

  if(graphics)  plots(density)

  return(density)
}


