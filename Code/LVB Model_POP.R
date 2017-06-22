# Bronwyn MacDonald
# May 9 2017

# Fitting von Bertalanffy Model to Chinook data
# For evaluation of size limits in ocean fisheries 

AL_DATA <- read.csv("CK Age-Length Data.csv")

# Code that installs packages if need be - taken from Forecast.R
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE)
    require(p, character.only = TRUE)
}

# Install packages
usePackage("R2jags")
usePackage("lattice")
usePackage("superdiag")
usePackage("FSA")
usePackage("dplyr")
usePackage("ggplot2")
usePackage("mcmcplots")
usePackage("nlstools")
usePackage("RColorBrewer")
# =================== Ludwig von Bertalanffy Length growth model ========================

Run_LVB <- function(Data){
  attach(Data)
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- Linf * (1-exp(-k * (Age[i] - t0)))
      Rep[i] ~ dnorm(eta[i], prec)
    }
    
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    sigma2 <- 1 /prec
    
    Linf ~ dunif(0, 10000)
    k     ~ dunif(0, 100)     
    t0    ~ dunif(-1000, 1000)
  }
  
   
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=1)
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=1)
  }

  data <- c("N", "Age", "Lmm")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  
  jags.LVB = jags.parallel(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)  	
  detach(Data)
  return(jags.LVB)

} # End Jags code


# =================== LOG ERROR ***** Ludwig von Bertalanffy Length growth model ========================

Run_LOG_LVB <- function(Data,Inits){
  Data$logLmm <- log(Data$Lmm)
  attach(Data)
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  N=length(Data$Lmm)
  
  LVBmodel <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- log(Linf * (1-exp(-k * (Age[i] - t0))))
     # Rep[i] ~ dnorm(eta[i], prec)
    }
    
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    sigma2 <- 1 /prec
    
    Linf ~ dunif(1, 10000)
    k     ~ dunif(0.001, 100)     
    t0    ~ dunif(-1000, 2)
  }
  
  
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=1)
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=1)
  }
#   
  data <- c("N", "Age", "logLmm")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  
  jags.LVB = jags.parallel(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)    
  detach(Data)
  return(jags.LVB)
  
} # End Jags code


# ========================== Normal LVB Model with stock index =======================

Run_Stock_LVB <- function(Data){
  nPop = length(unique(Data$Pop))
  attach(Data)
  N <- N[1]
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]])))
      Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0, 100) 
      t0[i]   ~ dunif(-1000, 1000)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(svTypical.=svTypical){
    list(
         Linf=rep(unlist(svTypical[[1]]), times=nPop), 
         k=rep(unlist(svTypical[[2]]),times=nPop), 
         t0=rep(unlist(svTypical[[3]]),times=nPop), 
         prec=rep(1,times=nPop)
    )
    list(
           Linf=rep(unlist(svTypical[[1]]), times=nPop), 
           k=rep(unlist(svTypical[[2]]),times=nPop), 
           t0=rep(unlist(svTypical[[3]]),times=nPop), 
           prec=rep(1,times=nPop)
         )
  }

  data <- c("N", "Age", "Lmm", "nPop", "POP_ID")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  return(jags.LVB)
  
} # End Jags code


# ========================== Normal LVB Model with stock index =======================

Run_Stock_LOG_LVB <- function(Data){
  nPop = length(unique(Data$Pop))
  Data$logLmm <- log(Data$Lmm)
  attach(Data)
  N <- N[1]
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- log(Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]]))))
      Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0.001, 100) 
      t0[i]   ~ dunif(-1000, 2)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]    
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(svTypical.=svTypical){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
  }
  
  data <- c("N", "Age", "logLmm", "nPop", "POP_ID")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  return(jags.LVB)
  
} # End Jags code

# ========================================================= Diagnostics ============================================ #

Diag.plots  <- function(mod,name,Data){
  mod.mcmc <- as.mcmc(mod)
  mod.mat <- as.matrix(mod.mcmc)  
  
  
  pdf(c(paste("Diag Plots",name,".pdf")))
  # ------------------ Prior Plots -----------------#
  
  if(POP==FALSE){
      par(mfrow=c(2,2))
      # Plot Linf priors
      x=mod.mat[,colnames(mod.mat)=="Linf"] ;   xx=as.matrix(x) ;
      plot( density(runif(10000000, 1, 10000)), main = c(expression(paste("Density of Linf"))), 
           xlab="", col="blue", lty="dashed")
      abline(v=xx, col="red") 
      # Plot k prior
      x=mod.mat[,colnames(mod.mat)=="k"] ;   xx=as.matrix(x) ;   
      plot(density(runif(10000000, 0.001, 100)), main = c(expression(paste("Density of k"))), 
           xlab="", col="blue", lty="dashed")
      abline(v=xx, col="red")  
      # Plot t prior
      x=mod.mat[,colnames(mod.mat)=="t0"] ;   xx=as.matrix(x) ; 
      plot(density(runif(10000000, -1000, 2)), main = c(expression(paste("Density of t0"))), 
           xlab="", col="blue", lty="dashed")
      abline(v=xx, col="red")  
      # Plot precision prior
      plot(density(rgamma(10000000, 0.0001, 0.0001)   , bw="nrd0" ), main = c(expression(paste("Density of sigma2"))), 
           xlab="", col="blue", lty="dashed")
  }
 
  # -------- Convergence Diagnostics --- #
 # summary(mod.mcmc)
  if(POP==TRUE){
      plot( density(runif(10000000, 1, 10000)), main = c(expression(paste("Density of Linf"))), 
            xlab="", col="blue", lty="dashed")
      plot(density(runif(10000000, 0.001, 100)), main = c(expression(paste("Density of k"))), 
           xlab="", col="blue", lty="dashed")
      plot(density(runif(10000000, -1000, 2)), main = c(expression(paste("Density of t0"))), 
           xlab="", col="blue", lty="dashed")
      plot(density(rgamma(10000000, 0.0001, 0.0001)   , bw="nrd0" ), main = c(expression(paste("Density of sigma2"))), 
           xlab="", col="blue", lty="dashed")
  }
  # Trace Plots
  if(ncol(mod.mat)<= 5)print(xyplot(mod.mcmc))
  if(ncol(mod.mat)>5){
    nplots <- ceiling(length(mod.mat[1,])/5)
    for(p in 1:(nplots-1)){
      print(xyplot(mod.mcmc[,((p*5)-4):(p*5)]))
    }
    xyplot(mod.mcmc[,(p*5+1):(length(mod.mcmc[[1]][1,])) ])
  } 
  print(densityplot(mod.mcmc, layout=c(2,3), aspect="fill", new=TRUE))
  print(gelman.plot(mod.mcmc))
  print(geweke.plot(mod.mcmc))
  superdiag(mod.mcmc, 250)
  wd <- getwd()
  mcmcplot(mod.mcmc, dir=wd)
  caterplot(mod.mcmc, parms = c("Linf", "k", "t0", "sigma2")) 
  autocorr.plot(mod.mcmc)
 dev.off()
}


# ================================================ Plot model Fit ================================================== #

# Prediction of Lengths
Predict <- function(Linf, k, t0, Age){
  Linf * (1-exp(-k * (Age - t0)))
}

Predict_POP <- function(Linf, k, t0, Age){
  Linf[POP_ID] * (1-exp(-k[POP_ID] * (Age - t0[POP_ID])))
}
  
Model_Plots <- function(mod, Plot_Data, name){
  pdf(c(paste("Fit Plot",name,".pdf")))
  cols <- brewer.pal(length(unique(Plot_Data$POP_ID)) ,"Dark2")
  
  mod.mcmc <- as.mcmc(mod)
  mod.mat <- as.matrix(mod.mcmc)  
  samp = sample(1:length(mod.mat[,1]), 1000)  # sample 1000 lines from parameters
  mod.sample = mod.mat[samp,]
  
  if(POP==FALSE){
    kpos = which(colnames(mod.sample)=="k")
    Lpos = which(colnames(mod.sample)=="Linf")
    tpos = which(colnames(mod.sample)=="t0")
  }
  if(POP==TRUE){
    kpos = which(colnames(mod.sample)=="k[1]")
    Lpos = which(colnames(mod.sample)=="Linf[1]")
    tpos = which(colnames(mod.sample)=="t0[1]")
  }
  pred.age <- seq(min(Plot_Data$Age),max(Plot_Data$Age), by=0.1)
  if(POP==FALSE){
    Pred <- sapply(pred.age, FUN=Predict, Linf= mod.sample[,Lpos], k=mod.sample[,kpos], t0=mod.sample[,tpos])
    Pred.plot <- data.frame(cbind(pred.age, 
                                  ave=apply(Pred, 2, FUN=quantile,probs=c(0.5)), 
                                  lower= apply(Pred, 2, FUN=quantile,probs=c(0.05)),
                                  upper= apply(Pred, 2, FUN=quantile,probs=c(0.95))))
    p=ggplot(aes(Age, Lmm), data=Plot_Data) + 
      geom_point() + 
      geom_line(color="blue",data=Pred.plot, aes(pred.age, ave )) +
      geom_line(color="red",data=Pred.plot, aes(pred.age, lower )) +
      geom_line(color="red",data=Pred.plot, aes(pred.age, upper ))
    print(p)
  }#end if POP==FALSE
  # For population models each element in the list Pred becomes predictions using the sampled paramters for a population
  if(POP==TRUE){
      Pred <- list()
      plot(Plot_Data$Age, Plot_Data$Lmm, xlab="Age", ylab="Length", cex=0.7)
      for(p in unique(Plot_Data$POP_ID)){
        Pred[[p]] <- sapply(pred.age, FUN=Predict, Linf= mod.sample[,(Lpos+p-1)], k=mod.sample[,(kpos+p-1)], t0=mod.sample[,(tpos+p-1)])
        lines(pred.age, apply(Pred[[p]], 2, FUN=quantile,probs=c(0.5)), col=cols[p])
        
        polygon(c(pred.age, rev(pred.age)), c(apply(Pred[[p]], 2, FUN=quantile,probs=c(0.5)), 
                                              rev(apply(Pred[[p]], 2, FUN=quantile,probs=c(0.05)))), col = alpha(cols[p],0.25), border = NA)
        polygon(c(pred.age, rev(pred.age)), c(apply(Pred[[p]], 2, FUN=quantile,probs=c(0.95)), 
                                              rev(apply(Pred[[p]], 2, FUN=quantile,probs=c(0.5)))), col = alpha(cols[p],0.25), border = NA)
       # lines(pred.age, apply(Pred[[p]], 2, FUN=quantile,probs=c(0.05)), lty="dashed", col=alpha(cols[p],0.5))
        #lines(pred.age, apply(Pred[[p]], 2, FUN=quantile,probs=c(0.95)), lty="dashed", col=alpha(cols[p],0.5))
      }
  }# end POP==TRUE 
  
 dev.off()
} # End Model.Plots




# ============================================ Select data to use ============================================= #


# Create data list

Get_CWT_Data <- function(AL_DATA){
  Data <- filter(AL_DATA, !is.na(CWT_Age))
  Data <- transform(Data, POP_ID=match(Pop, unique(Pop)))
  list(   
    N=length(Data$CWT_Age), 
    Age=Data$CWT_Age, 
    Lmm=Data$Est_FL,
    Pop=Data$Pop,
    POP_ID=Data$POP_ID
  )
}

Get_Scale_Data <- function(AL_DATA){
  Data <- filter(AL_DATA, !is.na(Scale_Agg_Age))
  Data <- transform(Data, POP_ID=match(Pop, unique(Pop)))
  list(   
    N=length(Data$Scale_Agg_Age), 
    Age=Data$Scale_Agg_Age, 
    Lmm=Data$Est_FL,
    Pop=Data$Pop,
    POP_ID=Data$POP_ID
  )
}

Get_Paired_CWT_Scale_Data <- function(AL_DATA){
  Data <- AL_DATA %>% filter(CWT_Age != "" & !is.na(Scale_Agg_Age))   # 6,924 samples
  Data <- transform(Data, POP_ID=match(Pop, unique(Pop)))
  list(   
    N=length(Data$Scale_Agg_Age),
    CWTAge=Data$CWT_Age,
    ScaleAge=Data$Scale_Agg_Age, 
    Lmm=Data$Est_FL,
    Pop=Data$Pop,
    POP_ID=Data$POP_ID
  )
}


# ======================  Run Both Freq and Bayes ==## ==========must run manually bc nls functions do not work inside functions!==== #


Run_Model <- function(Data=CWT_plot,Bayes_Data=CWT_Data, LOG=FALSE,POP=FALSE ,nboots=300){
  Data$logLmm <- log(Data$Lmm) # add log TL variable to crm
  Data <<- Data 
  
  if(POP==FALSE & LOG==FALSE) name = "Normal Model"
  if(POP==FALSE & LOG==TRUE)  name= "Log Model"
  if(POP==TRUE & LOG==FALSE) name = "Pop Normal Model"
  if(POP==TRUE & LOG==TRUE)  name= "Pop Log Model"
  
  pdf(c(paste("Freq Plots",name,".pdf")))
  
  if(POP==FALSE & LOG == FALSE){
    # Run frequentist estimation 
    svTypical <<- vbStarts( Lmm~Age, data = Data)
    vbTypical <<- Lmm ~ Linf*(1-exp(-K*(Age-t0)))
    fitTypical <- nls(vbTypical,data=Data,start=svTypical)
    bootTypical <- nlsBoot(fitTypical,niter=nboots) # niter should be nearer 1000
    print(confint(bootTypical,plot=TRUE))
    print(overview(fitTypical))
    
    # Fitted plot with confidence bounds showing uncertainty in the parameter estimates (blue) and residual uncertainty (red)
    ages2plot <- seq(min(Data$Age),max(Data$Age), by=0.1)
    fitPlot(fitTypical,xlab="Age",ylab="Total Length (mm)",cex=0.7,xlim=range(ages2plot),main="")
    LCI <- UCI <- LPI <- UPI <- numeric(length(ages2plot))
    ests <- bootTypical$coefboot
    for (i in 1:length(ages2plot)){ 
      pv <- ests[,"Linf"]*(1-exp(-ests[,"K"]*(ages2plot[i]-ests[,"t0"])))
      LCI[i] <- quantile(pv,0.025)
      UCI[i] <- quantile(pv,0.975)
      LPI[i] <- quantile(pv-bootTypical$rse,0.025)
      UPI[i] <- quantile(pv+bootTypical$rse,0.975)
    }
    lines(UCI~ages2plot,type="l",col="blue",lty=2)
    lines(LCI~ages2plot,type="l",col="blue",lty=2)
    lines(UPI~ages2plot,type="l",col="red",lty=2)
    lines(LPI~ages2plot,type="l",col="red",lty=2)      
    # Plot residuals
    print(residPlot(fitTypical))  
  } # end if LOG=FALSE
  
  if(POP==FALSE & LOG==TRUE){
    # LOG Model
    svTypicalM <<- vbStarts(Lmm~Age,data=Data,type="typical") # get starting values as usual
    vbTypical <<- vbFuns("typical") # get RHS of typical function
    fitTypicalM <- nls(logLmm~log(vbTypical(Age,Linf,K,t0)),data=Data,start=svTypicalM)   
    print(overview(fitTypicalM))
    
    # Cannot use bootstrap function with logged data o must do manually
    estsM <- data.frame(matrix(rep(NA,4), ncol=4,nrow=nboots)); colnames(estsM)=c("Linf","K","t0","rse")
    for(i in 1:nboots){
      samp <- sample(1:length(Data$Lmm), size=length(Data$Lmm), replace=T)
      sample <- (Data[samp,])
      fitsamp <- nls(logLmm~log(vbTypical(Age,Linf,K,t0)),data=sample,start=svTypicalM)
      estsM[i,1:3] <- coef(fitsamp)
      estsM[i,4] <- sqrt( sum(( fitted(fitsamp) - sample$logLmm)^2)/(length(sample$Lmm)-2) )
    }        
    # Parameter histogram plots
    np <- ncol(estsM[,1:3])
    rows <- round(sqrt(np))
    cols <- ceiling(sqrt(np))
    op <- graphics::par("mfrow")
    graphics::par(mfrow = c(rows, cols))
    conf.level=0.95
    cl <- function(x) stats::quantile(x, c((1 - conf.level)/2, 1 - (1 - conf.level)/2))
    res <- t(apply(estsM, 2, cl))
    for (i in 1:np) {
      h <- hist(~estsM[, i], xlab = colnames(estsM)[i])
      plotrix::plotCI(mean(estsM[, i]), y = 0.95 * max(h$counts), 
                      li = res[i, 1], ui = res[i, 2], err = "x", pch = 19, lwd=2, add = TRUE
      )
    }
    graphics::par(mfrow = op)
    # Fitted plot with confidence bounds showing uncertainty in the parameter estimates (blue) and residual uncertainty (red)
    ages2plot <- seq(min(Data$Age),max(Data$Age), by=0.1)
    fitPlot(fitTypicalM,xlab="Age",ylab="Log (Total Length (mm))",cex=0.7,xlim=range(ages2plot),main="")
    LCI <- UCI <- LPI <- UPI <- numeric(length(ages2plot))
    for (i in 1:length(ages2plot)){ 
      lnpv <- log(estsM[,"Linf"]*(1-exp(-estsM[,"K"]*(ages2plot[i]-estsM[,"t0"]))))
      LCI[i] <- quantile(lnpv,0.025)
      UCI[i] <- quantile(lnpv,0.975)
      LPI[i] <- quantile(lnpv-estsM$rse,0.025)
      UPI[i] <- quantile(lnpv+estsM$rse,0.975)
    }
    lines(UCI~ages2plot,type="l",col="blue",lty=2)
    lines(LCI~ages2plot,type="l",col="blue",lty=2)
    lines(UPI~ages2plot,type="l",col="red",lty=2)
    lines(LPI~ages2plot,type="l",col="red",lty=2)     
    # Plot residuals
    print(residPlot(fitTypicalM))     
  } # End if LOG = TRUE
  
  
  ### POPULATION BASED SECTION  
  if(POP==TRUE){
    if(LOG == FALSE){
      # Run frequentist estimation 
      svTypical <<- vbStarts( Lmm~Age, data = Data)
      if(POP==TRUE) svTypical <<- lapply(svTypical,rep,times=length(unique(Data$POP_ID)))
      vbTypical <<- Lmm ~ Linf*(1-exp(-K*(Age-t0)))
      if(POP==TRUE) vbTypical <<- Lmm ~ Linf[POP_ID]*(1-exp(-K[POP_ID]*(Age-t0[POP_ID])))
      fitTypical <- nls(vbTypical,data=Data,start=svTypical)
      print(overview(fitTypical))
      print(residPlot(fitTypical))
      #source('nlsBoot fix.R')
      Linf <- coef(fitTypical)[1:length(unique(Data$POP_ID))]
      k <- coef(fitTypical)[(length(unique(Data$POP_ID))+1) : (2*length(unique(Data$POP_ID)))] 
      t0 <- coef(fitTypical)[((2*length(unique(Data$POP_ID)))+1) : (3*length(unique(Data$POP_ID)))] 
    } # end if LOG=FALSE
    if(LOG==TRUE){
      # LOG Model
      svTypicalM <<- vbStarts(Lmm~Age,data=Data,type="typical") # get starting values as usual
      if(POP==TRUE) svTypicalM <<- lapply(svTypicalM,rep,times=length(unique(Data$POP_ID)))
      vbTypical <<- vbFuns("typical") # get RHS of typical function
      if(POP==TRUE) vbTypical <<- logLmm ~ log(Linf[POP_ID]*(1-exp(-K[POP_ID]*(Age-t0[POP_ID]))))
      fitTypicalM <- nls(vbTypical,data=Data,start=svTypicalM)   
      print(overview(fitTypicalM))
      print(residPlot(fitTypicalM))
      #source('nlsBoot fix.R')
      Linf <- coef(fitTypicalM)[1:length(unique(Data$POP_ID))]
      k <- coef(fitTypicalM)[(length(unique(Data$POP_ID))+1) : (2*length(unique(Data$POP_ID)))] 
      t0 <- coef(fitTypicalM)[((2*length(unique(Data$POP_ID)))+1) : (3*length(unique(Data$POP_ID)))] 
    }# if log = false
  
    # Plot model fit (cannot work bootstrap to get uncertainty)
    ages2plot <- seq(min(Data$Age),max(Data$Age), by=0.1)
    plot(Data$Age, Data$Lmm, ylab="Total Length (mm)",cex=0.7,xlim=range(ages2plot),main="", col=Data$POP_ID)
    for(i in 1:length(unique(Data$POP_ID))){
      lines(ages2plot, Predict(Linf[i], k[i], t0[i], Age=ages2plot), col=i)
    }
  }# end if POP==TRUE  
  dev.off() # turn off Freq Plots pdf  
  


# ---------------------------- END FREQUENTIST ---------------------------------- # 

# -------------------------------  Run BAYESIAN model -----------------------------------#
  
  if(POP==FALSE & LOG==FALSE) Model <- Run_LVB(Bayes_Data)
  if(POP==FALSE & LOG==TRUE)  LOG.model <- Run_LOG_LVB(Bayes_Data)
  if(POP==TRUE & LOG==FALSE)  Pop.model <- Run_Stock_LVB(Bayes_Data)
  if(POP==TRUE & LOG==TRUE)   Pop.LOG.model <- Run_Stock_LOG_LVB(Bayes_Data)  

    # Diagnostics
  if(POP==FALSE & LOG==FALSE) Diag.plots(Model,name=name, Data)
  if(POP==FALSE & LOG==TRUE)  Diag.plots(LOG.model,name=name, Data)
  if(POP==TRUE & LOG==FALSE)  Diag.plots(Pop.model, name=name, Data)
  if(POP==TRUE & LOG==TRUE)   Diag.plots(Pop.LOG.model, name=name, Data)
  # Model fit
  
  if(POP==FALSE & LOG==FALSE) Model_Plots(mod=Model, Plot_Data=Data, name=name)
  if(POP==FALSE & LOG==TRUE)  Model_Plots(mod=LOG.model, Plot_Data=Data, name=name)
  if(POP==TRUE & LOG==FALSE)  Model_Plots(mod=Pop.model, Plot_Data=Data, name=name)
  if(POP==TRUE & LOG==TRUE)   Model_Plots(mod=Pop.LOG.model, Plot_Data=Data, name=name)

  # Output
  if(POP==FALSE & LOG==FALSE) return(Model)
  if(POP==FALSE & LOG==TRUE)  return(LOG.model)
  if(POP==TRUE & LOG==FALSE)  return(Pop.model)
  if(POP==TRUE & LOG==TRUE)  return(Pop.LOG.model)
} # End Run_MOdel

# ===================================== END RUN_MODEL ========================================#


# ***************************  DATA ENTRY SECTION ****************************** #




# Paired CWT and Scale Data for comparison of LW relationship defined by eact (same data used for each analysis)
    Paired_Data <- Get_Paired_CWT_Scale_Data(AL_DATA) 
    Paired_plot <- as.data.frame(Paired_Data)
    op <- par(mfrow=c(1,3))
    hist(Paired_Data$Lmm~Paired_Data$POP_ID,col="grey"); hist(Paired_Data$CWTAge,col="grey");hist(Paired_Data$ScaleAge,col="grey")
    par(op)
    # Check error structure
    Summ_CWT <- Summarize(Lmm~CWTAge+POP_ID,data=Paired_plot,numdigs=1)
    print(Summ_CWT)
    plot( sd~mean, data=Summ_CWT, col=POP_ID,pch=19)
    # Check Scale age
    Summ_Scale <- Summarize(Lmm~ScaleAge+POP_ID,data=Paired_plot,numdigs=1)
    print(Summ_Scale)
    plot( sd~mean, data=Summ_Scale, col=POP_ID,pch=19)
    
    # Paired data with Age=CWT  age
    P.CWT_Data <- Paired_Data; P.CWT_Data$Age<-Paired_Data$CWTAge
    P.CWT_plot <- as.data.frame(P.CWT_Data)
    
    # Paired dat awith Age=Scale Age
    P.Scale_Data <- Paired_Data; P.Scale_Data$Age<-Paired_Data$ScaleAge
    P.Scale_plot <- as.data.frame(P.Scale_Data)

# CWT ONLY data
    CWT_Data <- Get_CWT_Data(AL_DATA)  # Initial Values 841.609230   1.160274   1.257773
    CWT_plot <- as.data.frame(CWT_Data)
    op <- par(mfrow=c(1,2))
    hist(CWT_Data$Lmm); hist(CWT_Data$Age);
    par(op)
    # Check error structure
    Summ <- Summarize(Lmm~Age,data=CWT_plot,numdigs=1)
    print(Summ)
    plot(Summ$mean, Summ$sd, pch=19)
    # ** If this shows increasing trend consider log transforming

# Scale ONLY Data
    Scale_Data <- Get_Scale_Data(AL_DATA) 
    Scale_plot <- as.data.frame(Scale_Data)
    op <- par(mfrow=c(1,2))
    hist(Scale_Data$Lmm); hist(Scale_Data$Age);
    par(op)
    # Check error structure
    Summ <- Summarize(Lmm~Age,data=Scale_plot,numdigs=1)
    print(Summ)
    plot(Summ$mean, Summ$sd, pch=19)
    
 
    
    
    
# =================================== RUNNING MODELS ============================ #
    
    
   # -------------------- RUN Paired CWT-based analysis ----------------------
    Paired.CWT.Model <- Run_Model(Data=P.CWT_plot, Bayes_Data=P.CWT_Data, LOG=FALSE, nboots=300, POP=FALSE)     # Normal Model
    Paired.CWT.Log.Model <- Run_Model(Data=P.CWT_plot, Bayes_Data=P.CWT_Data, LOG=TRUE, nboots=300, POP=FALSE)  # LOG Model

  # --------------------- RUN Paired Scale-based analysis ----------------------
    Paired.Scale.Model <- Run_Model(Data=P.Scale_plot, Bayes_Data=P.Scale_Data, LOG=FALSE, nboots=300, POP=FALSE)
    Paired.Scale.Log.Model <- Run_Model(Data=P.Scale_plot, Bayes_Data=P.Scale_Data, LOG=TRUE, nboots=300, POP=FALSE)

  # ------------------- RUN CWT ONLY Model ------------------
    CWT.Model <- Run_Model(Data=CWT_plot, Bayes_Data=CWT_Data, LOG=FALSE, nboots=300, POP=FALSE)
    CWT.Log.Model <- Run_Model(Data=CWT_plot, Bayes_Data=CWT_Data, LOG=TRUE, nboots=300, POP=FALSE)
    
  # ------------------ Run Scale age ONLY Model---------------  
    Scale.Model <- Run_Model(Data=Scale_plot, Bayes_Data=Scale_Data, LOG=FALSE, nboots=300, POP=FALSE)
    Scale.Log.Model <- Run_Model(Data=CWT_plot, Bayes_Data=CWT_Data, LOG=TRUE, nboots=300, POP=FALSE)

    
  # RUN POPULATION Model
   # Creating additional function called Run_Pop_Model which runs teh thing ut has different frequentist section
   Paired.POP.CWT.Model <- Run_Model(Data=P.CWT_plot, Bayes_Data=P.CWT_Data, LOG=FALSE, nboots=300, POP=TRUE)     # Normal Model  
   Paired.POP.Scale.Model <- Run_Model(Data=P.Scale_plot, Bayes_Data=P.Scale_Data, LOG=FALSE, nboots=300, POP=TRUE) # Normal Model   
    
   Paired.POP.LOG.CWT.Model <- Run_Model(Data=P.CWT_plot, Bayes_Data=P.CWT_Data, LOG=TRUE, nboots=300, POP=TRUE) 
   Paired.POP.LOG.Scale.Model <- Run_Model(Data=P.Scale_plot, Bayes_Data=P.Scale_Data, LOG=TRUE, nboots=300, POP=TRUE) 


# =============================== ANALYSIS and PLOTS =============================== #
    
# *********** Comparison parameter plots for CWT vs Scale ages
    # Normal Models
    CWT.mcmc <- as.mcmc(Paired.CWT.Log.Model)
    Scale.mcmc <- as.mcmc(Paired.Scale.Log.Model)
    
    op <- par(mfrow=c(2,2), mar=c(2,2,2,1))
    Linf_pos <- which(varnames(CWT.mcmc)=="Linf"); k_pos <- which(varnames(CWT.mcmc)=="k"); t0_pos <- which(varnames(CWT.mcmc)=="t0") 
    boxplot(unlist(CWT.mcmc[,Linf_pos]), unlist(Scale.mcmc[,Linf_pos]), main="Linf", names=c("CWT Data","Scale Data"), cex.label=0.7,cex.axis=0.7)
    boxplot(unlist(CWT.mcmc[,k_pos]), unlist(Scale.mcmc[,k_pos]), main="k", names=c("CWT Data","Scale Data"), cex.axis=0.7, cex.label=0.7)
    boxplot(unlist(CWT.mcmc[,t0_pos]), unlist(Scale.mcmc[,t0_pos]), main="t0", names=c("CWT Data","Stream Data"), cex.axis=0.7,cex.label=0.7)
    par(op)
    # ************************************************
    
# For POP Models
  # Do we need to bias correct scale data?
  # Parameter estiamtes are not sig different between CWT and scale data
  # *********** Comparison parameter plots for CWT vs Scale ages
 pdf("Comparison boxplots CWT vs. Scale POP LOG.pdf") 
 # Normal Models
  CWT.mat <- as.matrix(as.mcmc(Paired.POP.LOG.CWT.Model))
  Scale.mat <- as.matrix(as.mcmc(Paired.POP.LOG.Scale.Model))
  
  Linf_pos <- which(colnames(CWT.mat)=="Linf[1]"); k_pos <- which(colnames(CWT.mat)=="k[1]"); t0_pos <- which(colnames(CWT.mat)=="t0[1]") 
  
  L_plot <- rbind(melt(CWT.mat[,Linf_pos:(Linf_pos+length(unique(P.CWT_Data$POP_ID))-1) ]), 
                  melt(Scale.mat[,Linf_pos:(Linf_pos+length(unique(P.CWT_Data$POP_ID))-1) ]))
  L_plot$Parameter<-L_plot$Var2
  L_plot$Data <- c(rep("CWT", length(L_plot$Var1)/2), rep("Scale", length(L_plot$Var1)/2))
  ggplot(L_plot, aes(x = Parameter, y = value, fill = Data), xlab="") +
        geom_boxplot() +
        scale_fill_manual(values = c("yellow", "orange")) 
  
  k_plot <- rbind(melt(CWT.mat[,k_pos:(k_pos+length(unique(P.CWT_Data$POP_ID))-1) ]), 
                  melt(Scale.mat[,k_pos:(k_pos+length(unique(P.CWT_Data$POP_ID))-1) ]))  
  k_plot$Parameter<-k_plot$Var2
  k_plot$Data <- c(rep("CWT", length(k_plot$Var1)/2), rep("Scale", length(k_plot$Var1)/2))
  ggplot(k_plot, aes(x = Parameter, y = value, fill = Data), xlab="") +
    geom_boxplot() +
    scale_fill_manual(values = c("blue", "lightblue")) 

  t_plot <- rbind(melt(CWT.mat[,t0_pos:(t0_pos+length(unique(P.CWT_Data$POP_ID))-1) ]), 
                  melt(Scale.mat[,t0_pos:(t0_pos+length(unique(P.CWT_Data$POP_ID))-1) ]))  
  t_plot$Parameter<-t_plot$Var2
  t_plot$Data <- c(rep("CWT", length(t_plot$Var1)/2), rep("Scale", length(t_plot$Var1)/2))
  ggplot(t_plot, aes(x = Parameter, y = value, fill = Data), xlab="") +
    geom_boxplot() +
    scale_fill_manual(values = c("indianred4", "indianred1")) 

dev.off()
  # ************************************************



    # ******** Comparison model fit plots for paired CWT and Scale data
    CWT.mat <- as.matrix(CWT.mcmc)  
    Scale.mat <- as.matrix(Scale.mcmc)      
    samp = sample(1:length(CWT.mat[,1]), 1000)  # sample 1000 lines from parameters
    CWT.sample = CWT.mat[samp,]; Scale.sample = Scale.mat[samp,]
    kpos = which(colnames(CWT.sample)=="k")
    Lpos = which(colnames(CWT.sample)=="Linf")
    tpos = which(colnames(CWT.sample)=="t0")
    
    pred.age <- seq(min(Paired_Data$ScaleAge),max(Paired_Data$ScaleAge), by=0.1)
    Pred.CWT <- sapply(pred.age, FUN=Predict, Linf= CWT.sample[,Lpos], k=CWT.sample[,kpos], t0=CWT.sample[,tpos])
    Pred.Scale <- sapply(pred.age, FUN=Predict, Linf= Scale.sample[,Lpos], k=Scale.sample[,kpos], t0=Scale.sample[,tpos])
    Pred.plot <- data.frame(cbind(pred.age, 
                                  CWTave=apply(Pred.CWT, 2, FUN=quantile,probs=c(0.5)), 
                                  CWTlower= apply(Pred.CWT, 2, FUN=quantile,probs=c(0.05)),
                                  CWTupper= apply(Pred.CWT, 2, FUN=quantile,probs=c(0.95)),
                                  Scaleave=apply(Pred.Scale, 2, FUN=quantile,probs=c(0.5)), 
                                  Scalelower= apply(Pred.Scale, 2, FUN=quantile,probs=c(0.05)),
                                  Scaleupper= apply(Pred.Scale, 2, FUN=quantile,probs=c(0.95))            
                                  )
                            ) 
    op <- par(mar=c(3,3,1,1))
    plot(Lmm~ScaleAge, data=Paired_plot, axes=F,col="black", pch=19)
    axis(1, cex.axis=0.7); mtext(1,text="Age", cex=0.8, line=1.9)
    axis(2, cex.axis=0.7); mtext(2, text="Length", cex=0.8, line=2)
    points(Paired_plot$CWTAge,Paired_plot$Lmm, col="darkgrey", pch=19, cex=0.8)  
    lines(Pred.plot$pred.age, Pred.plot$Scaleave,  col="black",lwd=3)
    lines(Pred.plot$pred.age, Pred.plot$CWTave,  col="darkgrey",lwd=2)
 # **** End of Comparison model fit plot


 


