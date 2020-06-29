library('np')
options(np.messages=FALSE)
library('energy')
library('stats')
library('data.table')
library('boot')
# data generation mechanism used in the paper
# Returns a Matrix with n rows and columns (x, y, z_1, ...z_dim)
#x and y have to dimension real valued random variables
data.gen<-function(n,sigma_alt,dim){
	Data_to_write<- NULL
	z_prime = rnorm(n,0,sigma_alt)
	z = matrix(rnorm(n*dim,0,.3),nr=n, nc=dim)
	x = rnorm(n,z[,1],.2) + z_prime
	y = rnorm(n,z[,1],.2)  + z_prime
	Data_to_write=cbind(x,y,z)
	colnames(Data_to_write)[3:(dim+2)]= paste('z_', 1:dim,sep='')
	return(Data_to_write)
}
##############################################################################
##Calculating the Conditional distribution bandwidth objects and
## the \hatF.x_z(x_i|z_i), \hatF.y_z(y_i|z_i),\hatF.z(z_i)
##############################################################################

#F.x_z denotes F(X|Z) and Fbw.x_z denotes the  bandwidths used to estimate F.x_z

# F.z_hat is a matrix with columns  ( F(Z_d|Z_{1,..,(d-1)}, F(Z_{d-1}|Z_{1..(d-1)}), ..., F(Z_2|Z_1), F(Z_1))

Calc.F.cond<-function(x, y,z,z.dim, bwmet,bwobj=NULL){ # Conditional Distribution calculation at the data points
#bwobj should be of the kind that is returned from  Calc.F.cond. Use it in case you want to use previously calculated bandwdiths  (from other data points) to evaluate the estimate the conditional distribution functions from (x,y,z)

  if (bwmet=='fixed'){ # see np package for explanation of types of bw methods
    n<-length(x)
    Fbw.x_z<-rep(1.587*n^(-1/(4+1+z.dim)), z.dim+1)
    Fbw.y_z<-rep(1.587*n^(-1/(4+1+z.dim)),z.dim+1)
    temp <- Cond.dist.z(z,z.dim,bwmet)
  } else if (bwmet=='normal-reference'){
      Fbw.x_z <- npcdistbw(xdat=z, ydat=x,bwmethod=bwmet)
      Fbw.y_z <- npcdistbw(xdat=z, ydat=y,bwmethod=bwmet)
      temp <- Cond.dist.z(z,z.dim,bwmet)
  } else if (bwmet=='cv.ls'){
      if (is.null(bwobj)){  # on all this data x.dat is the Explanatory data and  y.dat is dependent data
        Fbw.x_z <- npcdistbw(xdat=z, ydat=x,bwmethod=bwmet)
        Fbw.y_z <- npcdistbw(xdat=z, ydat=y,bwmethod=bwmet)
        ## Calculating  F(Z_d|Z_{1,..,(d-1)} and F(Z_{d-1}|Z_{1..(d-1)})
        temp <- Cond.dist.z(z,z.dim,bwmet)
      }else { # if bwobj is not null we are using the previously estimated optimal bandwidths
        Fbw.x_z <-bwobj$Fbw.x_z
        Fbw.y_z <-bwobj$Fbw.y_z
        Fbw.z_z<- bwobj$Fbw.z_z
        ## Calculating  F(Z_d|Z_{1,..,(d-1)} and F(Z_{d-1}|Z_{1..(d-1)})
        temp <- Cond.dist.z(z,z.dim,bwmet, Fbw.z_z)
      }
  }
  ## Calculating  F(X|Z) and F(Y|Z)
  F.x_z <- fitted(npcdist(txdat=z,tydat=x,bws=Fbw.x_z)) ## using the data to train and fit
  F.y_z <- fitted(npcdist(txdat=z, tydat=y,bws=Fbw.y_z))
  return(list('F.x_z'=F.x_z,'F.y_z'=F.y_z,'Fbw.x_z'=Fbw.x_z,'Fbw.y_z'=Fbw.y_z,'Fbw.z_z'=temp$Fbw.z_z,'F.z_hat'=temp$F.z_hat))
}



Cond.dist.z<- function(z,z.dim, bwmet,bwobj=NULL){ #bwobj is vector of lists of type "condbandwidth" in the np package (the first element is of type "bandwidth" )
  n<-nrow(as.matrix(z))
  F.z_hat<-NULL# fitted values matrix (this is returned along with the bandwidth object used to get this fit)
    # on all this data x.dat is the Explanatory data and  y.dat is dependent data
  F.z_z<-vector("list",z.dim) # fitted values list (for temporary use)
  if (bwmet=='fixed'){
    Fbw.z_z<-vector("list",z.dim) # bandwidths
    Fbw.z_z[[1]]<-1.587*n^(-1/(2+z.dim))
     F.z_z[[1]]<-fitted(npudist(tdat=as.matrix(z)[,1], bws=Fbw.z_z[[1]]))
    if (z.dim >1){
      for (d in 2:z.dim) {
        z.ydat<-z[,d]#ydat<-Z_d
        z.xdat<-z[,1:(d-1)]  # xdat <-Z_{1,..,(d-1)}
        Fbw.z_z[[d]] <- rep(1.587*n^(-1/(4+d)), d)
        F.z_z[[d]] <-fitted(npcdist(txdat=z.xdat, tydat= z.ydat,bws=Fbw.z_z[[d]] )) ## same as adding exdat=z[,dim], eydat=z[,1:(dim-1)]
      }
    }
    for (temp.dim in 1:z.dim){
      F.z_hat<-cbind(F.z_hat,F.z_z[[temp.dim]])
    }
  }else if(bwmet=='normal-reference'){
     Fbw.z_z<-vector("list",z.dim) # bandwidths
      Fbw.z_z[[1]]<-npudistbw(dat=as.matrix(z)[,1], bwmethod='normal-reference')
       F.z_z[[1]]<-fitted(npudist(tdat=as.matrix(z)[,1], bws=Fbw.z_z[[1]]))
      if (z.dim >1){
        for (d in 2:z.dim) {
          z.ydat<-z[,d]#ydat<-Z_d
          z.xdat<-z[,1:(d-1)]  # xdat <-Z_{1,..,(d-1)}
          Fbw.z_z[[d]] <- npcdistbw(xdat=z.xdat, ydat=z.ydat, bwmethod='normal-reference')
          F.z_z[[d]] <- fitted(npcdist(bws=Fbw.z_z[[d]] )) ## same as adding exdat=z[,dim], eydat=z[,1:(dim-1)]
        }
      }
      for (temp.dim in 1:z.dim){
        F.z_hat<-cbind(F.z_hat,F.z_z[[temp.dim]])
      }
  }else if(bwmet=='cv.ls') {
    if (is.null(bwobj)) {
      Fbw.z_z<-vector("list",z.dim) # bandwidths
      Fbw.z_z[[1]]<-npudistbw(dat=as.matrix(z)[,1])
       F.z_z[[1]]<-fitted(npudist(tdat=as.matrix(z)[,1], bws=Fbw.z_z[[1]]))
      if (z.dim >1){
        for (d in 2:z.dim) {
          z.ydat<-z[,d]#ydat<-Z_d
          z.xdat<-z[,1:(d-1)]  # xdat <-Z_{1,..,(d-1)}
          Fbw.z_z[[d]] <- npcdistbw(xdat=z.xdat, ydat=z.ydat, bwmethod='cv.ls')
          F.z_z[[d]] <- fitted(npcdist(bws=Fbw.z_z[[d]] )) ## same as adding exdat=z[,dim], eydat=z[,1:(dim-1)]
        }
      }
      for (temp.dim in 1:z.dim){
        F.z_hat<-cbind(F.z_hat,F.z_z[[temp.dim]])
      }
    }else{
      Fbw.z_z<-bwobj
      F.z_z[[1]]<-fitted(npudist(tdat=as.matrix(z)[,1], bws=Fbw.z_z[[1]]))
      if (z.dim >1){
        for (d in 2:z.dim) {
          z.ydat<-z[,d]#ydat<-Z_d
          z.xdat<-z[,1:(d-1)]  # xdat <-Z_{1,..,(d-1)}
          F.z_z[[d]] <- fitted(npcdist(txdat=z.xdat, tydat= z.ydat,bws=Fbw.z_z[[d]] )) # using the given bandwidth and reestimating F.z at the new values
        }
      }
      for (temp.dim in 1:z.dim){
        F.z_hat<-cbind(F.z_hat,F.z_z[[temp.dim]])
      }
    }
  }
    return(list('Fbw.z_z'=Fbw.z_z, 'F.z_hat'=F.z_hat))
}





##############################################################################
# Test-statistic on multivariate STANDARD normality with optional Freq vector
##############################################################################

std.mvnorm.e <-function(x, freq.rep=NULL){
    # E-statistic for multivariate normality
    if (is.vector(x) & is.null(freq.rep)==F) return(std.normal.e(x,freq.rep))
    if (is.vector(x) & is.null(freq.rep)==T) return(std.normal.e(x))
    if (is.null(freq.rep)==F) x<-as.matrix(x)[rep(1:nrow(as.matrix(x)), freq.rep),]
    n <- nrow(x)
    d <- ncol(x)
    if (n < 2) return(std.normal.e(x))
    y <- x
    if (any(!is.finite(y))) return (NA)
    stat <- 0
    e <- .C("mvnEstat", y = as.double(t(y)), byrow = as.integer(TRUE),
            nobs = as.integer(n), dim = as.integer(d),
            stat = as.double(stat), PACKAGE = "energy")$stat
    e
}

std.normal.e <-function(x,freq.rep=NULL){
  if (is.null(freq.rep)==F) x<-x[rep(1:length(x), freq.rep)]
  # E-statistic for multivariate normality
   x <- as.vector(x)
   y <- sort(x)
   n <- length(y)
   if (y[1] == y[n]) return (NA)
   K <- seq(1 - n, n - 1, 2)
   e <- 2 * (sum(2 * y * pnorm(y) + 2 * dnorm(y)) - n / sqrt(pi) - mean(K * y))
   e
}
emp.mvnorm.e<-function(x, freq.rep=NULL){
  if (is.null(freq.rep)==F) x<-as.matrix(x)[rep(1:nrow(as.matrix(x)), freq.rep),]
  return (mvnorm.e(x))
}
emp.normal.e<-function(x, freq.rep=NULL){
  x<-as.vector(x)
  if (is.null(freq.rep)==F) x<-x[rep(1:length(x), freq.rep)]
  return (normal.e(x))
}

# ##############################################################################
# ALL CV
##############################################################################



doing.boot.F.hat.emp.z.cv<-function(Dist.calc,bwmet,x,y,z,n ,z.dim,accuracy,R = 999){
  stat.boot<-rep(0,R)
  # creating the boot sample for z by sampling from its empirical
  z<-as.matrix(z)
  index.z<-sample(1:dim(z)[1],size=R*n, replace = T)
  samp.z=z[index.z, ]
  samp.z=as.matrix(samp.z)
  # evaluating F.hat.x_z and F.hat.x_z at grid points

  zdata<-z[rep(1:n, each = accuracy),]
  rownames(zdata)<-NULL
  xvalues<-seq(min(x), max(x), length.out=accuracy)
  yvalues<-seq(min(y), max(y), length.out=accuracy)

  eval.F.x_z <- fitted(npcdist(exdat=zdata,  eydat= rep(xvalues,times=n), txdat=z, tydat=x, bws=Dist.calc$Fbw.x_z)) ## is this right?? yes
  eval.F.y_z <- fitted(npcdist(exdat=zdata,  eydat= rep(yvalues,times=n), txdat=z, tydat=y, bws=Dist.calc$Fbw.y_z))

  eval.mat.F.x_z<-matrix(eval.F.x_z,nrow=n, ncol=accuracy,byrow=T)
  eval.mat.F.y_z<-matrix(eval.F.y_z,nrow=n, ncol=accuracy,byrow=T)
  #sampling boot.x and boot.y from F.inv.hat.x_z and F.inv.hat.y_z
  unif.data<- matrix(runif(n*R*2,0,1),nrow=n*R, ncol=2)
  boot.x<-rep(0,n*R )
  boot.y<- rep(0,n*R )
  for (i in 1:n){
    temp.unif<-unif.data[index.z==i,]
    if (is.vector(temp.unif)) temp.unif<-matrix(temp.unif,ncol=2)
    boot.x[index.z==i]=approx(eval.mat.F.x_z[i,],xvalues,xout=temp.unif[,1], method='linear',rule=2)$y
    boot.y[index.z==i]=approx(eval.mat.F.y_z[i,],yvalues,xout=temp.unif[,2], method='linear',rule=2)$y
  }
  print(paste(R, 'bootstrap samples obtained'))

  # using the boot data to evaluate the distribution of the statistics
  for (boot.no in 1:R){
    ind.boot <- (boot.no-1)*n +1:n
    if(boot.no%%25==0) print(paste('At bootstrap iteration', boot.no, 'of' , R,sep=' '))
    temp.boot<-Calc.F.cond(boot.x[ind.boot],boot.y[ind.boot],samp.z[ind.boot,],z.dim,'cv.ls')
    F.hat.temp<-cbind(as.matrix(temp.boot$F.x_z),as.matrix(temp.boot$F.y_z),as.matrix(temp.boot$F.z_hat))
    F.hat.temp[F.hat.temp==1]<-1-10^(-3)
    F.hat.temp[F.hat.temp==0]<-10^(-3)
    stat.boot[boot.no]<-std.mvnorm.e(qnorm(F.hat.temp))
  }
  return(stat.boot)
}

boot.F.hat.emp.z.cv<-function(Dist.calc,x,y,z,n ,z.dim,accuracy,R = 999){
    print(paste('Starting bootstrap'))
    d=2+ncol(z)
    stat.boot<-doing.boot.F.hat.emp.z.cv(Dist.calc,'cv.ls',x,y,z,n ,z.dim,accuracy,R)
    temp<-cbind(Dist.calc$F.x_z, Dist.calc$F.y_z, Dist.calc$F.z_hat)
    colnames(temp)=c('F.x_z', 'F.y_z', paste('F.hat.z',1:z.dim,sep=''))
    temp[temp==1]<-1-10^(-3)
    temp[temp==0]<-10^(-3)
    stats<-std.mvnorm.e(qnorm(temp))
    p <- 1 - mean(stat.boot < stats)
  e <- list(statistic = stats,
              p.value = p,
              method = "Cond Indep test: p-values by inverting F_hat to get bootstrap samples",
              bandwidth.method = 'least-squares cross-validation, see "np" package',
              data.desrip = paste("dimension of Z is ", z.dim,", sample size ", n, ", dimension of (X,Y,Z) is ", d, ", bootstrap replication number", R, sep = ""),
              bootstrap.stat.values=stat.boot)
    e
}



npresid.statistics<-function(data.to.work,z.dim){
  data.to.write.list<-data.to.work
  ##############################################################################
  x<- data.to.work[,1]
  data.to.work<-data.to.work[,-1]
  y<-data.to.work[,1]
  data.to.work<-data.to.work[,-1]
  z<-as.matrix(data.to.work)[,1:z.dim]
  Dist.calc<-Calc.F.cond(x, y, z, z.dim,bwmet='cv.ls')
  temp<-cbind(Dist.calc$F.x_z, Dist.calc$F.y_z, Dist.calc$F.z_hat)
  colnames(temp)=c('F.x_z', 'F.y_z', paste('F.hat.z',1:z.dim,sep=''))
  temp[temp==1]<-1-10^(-3)
  temp[temp==0]<-10^(-3)
  stats<-std.mvnorm.e(qnorm(temp))
  return(stats)
}

npresid.boot<-function(data.to.work,z.dim,accuracy=500,boot.replic=999){
	data.to.write.list<-data.to.work
	##############################################################################
	x<- data.to.work[,1]
	data.to.work<-data.to.work[,-1]
	y<-data.to.work[,1]
	data.to.work<-data.to.work[,-1]
	z<-as.matrix(data.to.work)[,1:z.dim]
	Dist.calc<-Calc.F.cond(x, y, z, z.dim,bwmet='cv.ls')
	boot.data <-boot.F.hat.emp.z.cv(Dist.calc,x,y,z,n ,z.dim,accuracy,boot.replic)
	new.data<-c(boot.data, list(cond.dist.obj=Dist.calc))
	return(new.data)
}