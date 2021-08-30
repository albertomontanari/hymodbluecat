################################################################
### Hymod calibration
################################################################

hymod.par=function(param.ini,area,tdelta,e,p,nstep=length(p),qoss,qinitial=0,lower=c(10,0.1,0,1,1),upper=c(400,10,0.9,1000,1000),itermax=100,control=list(factr=1e1,fnscale=0.01,parscale=c(100,1,1,1,1)),opt="DEoptim",plot=T)
{
#Choose between DEoptim and optim
if (opt=="DEoptim")
  {cat("Launching optimisation with R function DEoptim - please wait for calibration to be completed", fill=T)
  pr1=DEoptim(hymod.eff,lower=lower,upper=upper,control = DEoptim.control(itermax=itermax),args=list(area=area,tdelta=tdelta,e=e,p=p,qoss=qoss,qinitial=qinitial,nstep=nstep))} else
  {cat("Launching optimisation with R function optim - please wait for calibration to be completed", fill=T)
  pr1=optim(param.ini,hymod.eff,method="L-BFGS-B",lower=lower,upper=upper,control=control,args=list(area=area,tdelta=tdelta,e=e,p=p,qoss=qoss,qinitial=qinitial,nstep=nstep))}
if (opt=="DEoptim")
  bpar=pr1$optim$bestmem else bpar=pr1$par
  qsim=hymod.sim(bpar,area,tdelta,e,p,qinitial=qinitial)
  pr1$qsim=qsim$q_tot
  pr1$qoss=qoss
  #Compute efficiency
  eff=1-sum((pr1$qsim-pr1$qoss)^2)/sum((pr1$qoss-mean(pr1$qoss))^2)
  pr1$eff=eff
  #Make scatterplot
  if (plot==T)
    {
    plot(pr1$qsim,pr1$qoss,xlim=c(0,max(cbind(pr1$qsim,pr1$qoss))),ylim=c(0,max(cbind(pr1$qsim,pr1$qoss))),xlab="Simulated data",ylab="Observed data",col="red")
    grid()
    abline(0,1,lwd=2)
    legend("topleft",legend=bquote(Efficiency== .(signif(eff,digit=2))),cex=1.3)
    }
  pr1$par=bpar
  return(pr1)
}


################################################################
### Hymod to compute negative efficiency: FORTRAN 95
################################################################

hymod.eff<-function(param,args)
{
#Extract arguments from the argument list
area=args$area
tdelta=args$tdelta
e=args$e
p=args$p
qoss=args$qoss
qinitial=args$qinitial
nstep=args$nstep
#Define qout
qout<-rep(999,nstep)
#Define eff
eff=0
#Call the fortran subroutine to run Hymod
qsimf<-.Fortran("hymodfortran", 
                    param=as.double(param),
                    area=as.double(area),
                    tdelta=as.integer(tdelta),
                    e=as.double(e),
                    p=as.double(p),
                    nstep=as.integer(nstep),
                    qinitial=as.double(qinitial),
                    evapt=as.double(p),
                    qtslow=as.double(p),
                    qtquick=as.double(p),
                    qtot=as.double(p),
                    PACKAGE="hymodbo")
#Compute the efficiency                    
qsimf$qtot[qsimf$qtot>(max(qoss)*2)]=max(qoss)*2
eff=-1+sum((qsimf$qtot-qoss)^2)/sum((qoss-mean(qoss))^2)
return(eff)
}

################################################################
### Hymod to return simulated flow with uncertainty: FORTRAN 95
################################################################

hymod.sim<-function(param,area,tdelta,e,p,resultcalib=NULL,nstep=length(p),qinitial=0,bluecat=F,predsmodel="avg",empquant=F,siglev=0.2,m=100,m1=80,paramd=c(0.1,1,10,NA),lowparamd=c(0.001,0.001,0.001,0),upparamd=c(1,5,20,NA),NSeff=F,qoss=NULL,plot=F,dataplot="dmodel",cpptresh=0)
{
#Define qout
qout<-rep(999,nstep)
#Call the fortran routine to run Hymod
qsimf<-.Fortran("hymodfortran", 
                    param=as.double(param),
                    area=as.double(area),
                    tdelta=as.integer(tdelta),
                    e=as.double(e),
                    p=as.double(p),
                    nstep=as.integer(nstep),
                    qinitial=as.double(qinitial),
                    evapt=as.double(p),
                    qtslow=as.double(p),
                    qtquick=as.double(p),
                    qtot=as.double(p),
                    PACKAGE="hymodbluecat")
#Choose whether to run Bluecat uncertainty estimation. If yes, argument resultcalib must be provided
if (bluecat==T)
  {
  if(is.null(resultcalib)==T) 
    {
    cat("Error: argument resultcalib is NULL. Exiting",fill=T)
    return()
    }
  #Variable zeta serves for the diagnostic plots if required
  zeta=rep(0,nstep)
  #Variable aux is used to order the simulated data in ascending order and later to put back stochastic predictions in chronological order
  aux=sort(qsimf$qtot,index.return=T)$ix
  #sortsim contains the simulated data in ascending order
  sortsim=sort(qsimf$qtot)
  #nstep1 is the length of the calibration data set
  nstep1=length(resultcalib$qsim)
  #Variable aux2 is used to order the simulated calibration data in ascending order and to order observed calibration data according to ascending order of simulated calibration data
  aux2=sort(resultcalib$qsim,index.return=T)$ix
  #Ordering simulated calibration data in ascending order
  sortcalibsim=sort(resultcalib$qsim)
  #Ordering observed calibration data in ascending order of simulated calibration data
  qossc=resultcalib$qoss[aux2]
  
  #Find the vector of minimal quantities as computed below. It serves to identify the range of observed data to fit
  vectmin=pmin(rep(0.5,nstep1)+(nstep1-seq(1,nstep1))*0.5/m/2,rep(1,nstep1))
  #Find the vector of minimal quantities as computed below. It serves to identify the range of observed data to fit
  vectmin1=floor(pmin(rep(m,nstep1),(seq(1,nstep1)-1),rev((seq(1,nstep1)-1))/vectmin))

  #Defines the vectors of stochastic prediction and confidence bands
  medpred=rep(0,nstep)
  infpred=rep(0,nstep)
  suppred=rep(0,nstep)

  #Definition of auxiliary variable icount, to be used only to print on screen the progress of computation 
  icount=1

  cat("Bluecat uncertainty estimation. This will take some time",fill=T)
  cat("Computing the mean stochastic prediction",fill=T)
  for (i in 1:nstep)
    {
    #Routine to show the progress of the computation
    if (i/nstep*100 > 10*icount) 
      {
      cat(paste(icount*10),"% ",sep="")
      icount=icount+1
      }
    if (i/nstep*100 > 99 && icount==10)
      {
      cat("100%",fill=T)
      icount=icount+1
      }
    #Finding the simulated data in the calibration period closest to the data simulated here. Values are multiplied by one million to avoid small differences leading to multiple matches
    indatasimcal=Closest(sortcalibsim*10^6,sortsim[i]*10^6, which = T, na.rm = T)
    #Second workaround to avoid multiple matches
    indatasimcal1=indatasimcal[length(indatasimcal)]
    #Define the end of the range of the observed data to fit
    aux1=indatasimcal1-vectmin1[indatasimcal1]+(1+vectmin[indatasimcal1])*vectmin1[indatasimcal1]
    #Puts a limit to upper index of conditioned vector of observed data
#    cat(i,indatasimcal1,vectmin1[indatasimcal1],vectmin[indatasimcal1],aux1,nstep,fill=T)
    if(aux1>nstep1) aux1=nstep1
    #Define the start of the range of the observed data to fit
    aux2=indatasimcal1-vectmin1[indatasimcal1]
    #Puts a limit to lower index of conditioned vector of observed data
    if(aux2<1) aux2=1 
    #Compute mean stochastic prediction
    if(predsmodel=="avg") medpred[i]=mean(qossc[aux2:aux1])
    if(predsmodel=="mdn") medpred[i]=median(qossc[aux2:aux1])
    }
  #Put back medpred in chronological order
  medpred=medpred[order(aux)]
  #Choose between empirical quantiles and k-moments quantiles
  if(empquant==F)
    {
    #Estimation of ph and pl - orders of the k-moments for upper and lower tail
    #Fitting of PBF distribution on the sample of mean stochastic prediction
    #Update the initial values of the PBF parameter xl
    paramd[4]=0.5*min(medpred)
    upparamd[4]=0.9*min(medpred)
    #Definition of the number of k-moments to estimate on the sample of mean stochastic prediction to fit the PBF distribution
    m2=seq(0,m1)
    #Definition of the order p of the k-moments to estimate on the sample of mean stochastic prediction to fit the PBF distribution
    ptot=nstep^(m2/m1)
    #Estimation of k-moments for each order. k-moments are in kp and kptail
    Fxarr1=rep(0,nstep)
    kp=rep(0,(m1+1))
    kptail=rep(0,(m1+1))
    for(ii in 1:(m1+1))
      {
      p1=ptot[ii]
      for(iii in 1:nstep)
        {
        if(iii<p1) c1=0 else if(iii<p1+1 || abs(c1)<1e-30)
          {
          c1=exp(lgamma(nstep-p1+1)-lgamma(nstep)+lgamma(iii)-lgamma(iii-p1+1)+log(p1)-log(nstep))} else 
          c1=c1*(iii-1)/(iii-p1)
          Fxarr1[iii]=c1
        }
        kp[ii]=sum(sort(medpred)*Fxarr1)
        kptail[ii]=sum(rev(sort(medpred))*Fxarr1)
      }
      #End estimation of k-moments
      #Fitting of PBF distribution by using DEoptim with default parameter bounds. Make sure
      #DEoptim is installed and loaded in the library. If it fails change the default values for
      #lowparamd and upparamd
      PBFparam=DEoptim(fitPBF,lower=lowparamd,upper=upparamd,ptot=ptot,kp=kp,kptail=kptail)
      paramd=unname(PBFparam$optim$bestmem)
      #Recomputation of lambda1 and lambdainf with the calibrated parameters
      lambda1k=(+1+(beta(1/paramd[2]/paramd[1]-1/paramd[2],1/paramd[2])/paramd[2])^paramd[2])^(1/paramd[2]/paramd[1])
      lambda1t=1/(1-(+1+(beta(1/paramd[2]/paramd[1]-1/paramd[2],1/paramd[2])/paramd[2])^paramd[2])^(-1/paramd[2]/paramd[1]))
      lambdainfk=gamma(1-paramd[1])^(1/paramd[1])
      lambdainft=gamma(1+1/paramd[2])^(-paramd[2])
      #Computation of orders ph and pl of -kmoments
      ph=1/(lambdainfk*siglev/2)+1-lambda1k/lambdainfk
      pl=1/(lambdainft*siglev/2)+1-lambda1t/lambdainft
    }
  #Definition of auxiliary variable icount, to be used only to print on screen the progress of computation 
  icount=1
  #Routine for computing upper and lower confidence bands
  cat("Computing prediction confidence bands",fill=T)
  for (i in 1:nstep)
    {
    #Finding the simulated data in the calibration period closest to the data simulated here. Values are multiplied by one million to avoid small differences leading to multiple matches
    indatasimcal=Closest(sortcalibsim*10^6,sortsim[i]*10^6, which = T, na.rm = T)
    #Second workaround to avoid multiple matches
    indatasimcal1=indatasimcal[length(indatasimcal)]
    #Define the end of the range of the observed data to fit
    aux1=indatasimcal1-vectmin1[indatasimcal1]+(1+vectmin[indatasimcal1])*vectmin1[indatasimcal1]
    #Puts a limit to upper index of conditioned vector of observed data
    if(aux1>nstep1) aux1=nstep1
    #Define the start of the range of the observed data to fit
    aux2=indatasimcal1-vectmin1[indatasimcal1]
    #Puts a limit to lower index of conditioned vector of observed data
    if(aux2<1) aux2=1 
    #Defines the size of the data sample to fit
    count=aux1-aux2+1
    if(empquant==F)
      {
      #Routine to show the progress of the computation
      if (i/nstep*100 > 10*icount) 
        {
        cat(paste(icount*10),"% ",sep="")
        icount=icount+1
        }
      if (i/nstep*100 > 99 && icount==10)
        {
        cat("100%",fill=T)
        icount=icount+1
        }
      #Estimation of k moments for each observed data sample depending on ph (upper tail) and pl (lower tail)
      Fxpow1=ph-1
      Fxpow2=pl-1
      Fxarr1=rep(0,count)
      Fxarr2=rep(0,count)
      for(ii in 1:count)
        {
        if(ii<ph) c1=0 else if(ii<ph+1 || abs(c1)<1e-30)
          {
          c1=exp(lgamma(count-ph+1)-lgamma(count)+lgamma(ii)-lgamma(ii-ph+1)+log(ph)-log(count))} else 
          c1=c1*(ii-1)/(ii-ph)
        if(ii<pl) c2=0 else if(ii<pl+1 || abs(c2)<1e-30)
          {
          c2=exp(lgamma(count-pl+1)-lgamma(count)+lgamma(ii)-lgamma(ii-pl+1)+log(pl)-log(count))} else 
          c2=c2*(ii-1)/(ii-pl)
          Fxarr1[ii]=c1
          Fxarr2[ii]=c2
        }
      #End of k-moments estimation
      #Computation of confidence bands
      suppred[i]=sum(sort(qossc[aux2:aux1])*Fxarr1)
      infpred[i]=sum(rev(sort(qossc[aux2:aux1]))*Fxarr2)
      #Do not compute confidence bands with less than 3 data points
      if(count<3) 
      {
      suppred[i]=NA
      infpred[i]=NA
      }
      } else
      #Empirical quantile estimation - This is much easier
      {
      #Routine to show the progress of the computation
      if (i/nstep*100 > 10*icount) 
        {
        cat(paste(icount*10),"% ",sep="")
        icount=icount+1
        }
      if (i/nstep*100 > 99 && icount==10)
        {
        cat("100%",fill=T)
        icount=icount+1
        }
      #Computation of the position (index) of the quantiles in the sample of size count
      eindexh=ceiling(count*(1-siglev/2))
      eindexl=ceiling(count*siglev/2)
      #Computation of confidence bands
      suppred[i]=sort(qossc[aux2:aux1])[eindexh]
      infpred[i]=sort(qossc[aux2:aux1])[eindexl]
      if(count<3) 
      {
      suppred[i]=NA
      infpred[i]=NA
      }
      }
      #If plot=T compute the data to draw the diagnostic plots
      if(plot==T && nstep1>20 && is.null(qoss)==F)
        {
        #Defines the z vector
        qossaux=qoss[aux]
        #Finds the index in sorted vector of data defining the conditional distribution that is closest to the considered observed value
        indataosssim=Closest(sort(qossc[aux2:aux1])*10^6,qossaux[i]*10^6, which = T, na.rm = T)
        #Workaround to avoid multiple matches
        indataosssim1=indataosssim[length(indataosssim)]
        #Finds probability with Weibull plotting position, only if there are at least three data points
        if(count>=3) zeta[i]=indataosssim1/(count+1)
        }
    }
  #Put confidence bands back in chronological order
  infpred=infpred[order(aux)]
  suppred=suppred[order(aux)]
  }
#Diagnostic activated only if bluecat=T, plot=T, nstep1>20 and qoss is different from NA
if(bluecat==T && plot==T && nstep1<=20) cat("Cannot perform diagnostic if simulation is shorter than 20 time steps",fill=T)
if(bluecat==T && plot==T && is.null(qoss)==T) cat("Cannot perform diagnostic if observed flow is not available",fill=T)
if(bluecat==T && plot==T && nstep1>20 && is.null(qoss)==F)
  {
  #Plotting the diagnostic plots and scatterplot
  #Datapoints with observed flow lower than cpptresh are removed
  #Preparing data for D-model plot
  sortdata=sort(resultcalib$qoss)
  sortdata=sort(qsimf$qtot[qoss>cpptresh])
  sortdata1=sort(medpred[qoss>cpptresh])
  zQ=rep(0,101)
  z=seq(0,1,0.01)
  fQ=sortdata[ceiling(z[2:101]*length(sortdata))]
  fQ=c(0,fQ)
  zQ[1]=0
  #Preparing data for S-model plot
  zQ1=rep(0,101)
  z1=seq(0,1,0.01)
  fQ1=sortdata1[ceiling(z1[2:101]*length(sortdata1))]
  fQ1=c(0,fQ1)
  zQ1[1]=0
  qoss2=qoss[qoss>cpptresh]
  for(i in 2:101)
    {
    zQ[i]=(length(qoss2[qoss2<fQ[i]]))/length(sortdata)
    zQ1[i]=(length(qoss2[qoss2<fQ1[i]]))/length(sortdata1)
    }
  #Combining 4 plots in the same window
  if(is.null(dev.list())==F) dev.off()
  dev.new(width=12,height=10)
  par(mar=c(5,6,3,3))
  par(mfrow=c(2,2))  
  # First plot
  plot(z,zQ,ylim=c(0,1),xlim=c(0,1),type="b",xlab="z",ylab=expression('F'[q]*'(F'[Q]*''^-1*'(z))')) 
  points(z1,zQ1,ylim=c(0,1),xlim=c(0,1),pch=19,col="red") 
  legend(0,1,inset=0,legend=c("D-model", "S-model"),col=c("black", "red"),cex=1.2,pch=c(1,19))
  abline(0,1)  
  grid()
  title("Predictive cumulative distribution plot")
  # Second plot
  plot(sort(zeta[zeta!=0]),ppoints(zeta[zeta!=0]),ylim=c(0,1),xlim=c(0,1),type="b",xlab="z",ylab=expression('F'[z]*'(z)')) 
  abline(0,1)  
  grid()
  title("Combined probability-probability plot")
  aux4=sort(medpred,index.return=T)$ix
  #sortmedpred contains the data simulated by stochastic model in ascending order
  sortmedpred=sort(medpred)
  #Ordering observed calibration data and confidence bands in ascending order of simulated data by the stochastic model
  sortmedpredoss=qoss[aux4]
  sortsuppred=suppred[aux4]
  sortinfpred=infpred[aux4]
  # Third plot
  plot(sortmedpred,sortmedpredoss,xlim=c(0,max(cbind(sortmedpred,sortmedpredoss))),ylim=c(0,max(cbind(sortmedpred,sortmedpredoss))),xlab="Data simulated by stochastic model",ylab="Observed data",col="red")
  #Adding confidence bands. Not plotted for higher values of flow to avoid discontinuity
  aux5=nstep-10
  lines(sortmedpred[1:aux5],sortsuppred[1:aux5],col="blue")
  lines(sortmedpred[1:aux5],sortinfpred[1:aux5],col="blue")
  abline(0,1)
  grid()
  #Compute the efficiency of the simulation and put it in the plots
  eff=1-sum((medpred-qoss)^2)/sum((qoss-mean(qoss))^2)
  eff1=1-sum((qsimf$qtot-qoss)^2)/sum((qoss-mean(qoss))^2)
  legend("topleft",inset=0,legend=bquote("Efficiency S-model ="~.(signif(eff,digit=2))),cex=1.3)
  title("Scatterplot S-model predicted versus observed data")
  # Fourth plot
  plot(qsimf$qtot[aux4],sortmedpredoss,xlim=c(0,max(cbind(qsimf$qtot[aux4],sortmedpredoss))),ylim=c(0,max(cbind(qsimf$qtot[aux4],sortmedpredoss))),xlab="Data simulated by deterministic model",ylab="Observed data",col="red")
  abline(0,1)
  grid()
  #Compute the efficiency of the simulation and put it in the plots
  eff=1-sum((medpred-qoss)^2)/sum((qoss-mean(qoss))^2)
  eff1=1-sum((qsimf$qtot-qoss)^2)/sum((qoss-mean(qoss))^2)
  legend("topleft",inset=0,legend=bquote("Efficiency D-model ="~.(signif(eff1,digit=2))),cex=1.3)
  title("Scatterplot D-model predicted versus observed data")

  
  #Computing percentages of points lying outside the confidence bands by focusing on observed data greater than cpptresh
  cat("Percentage of points lying above the upper confidence limit=",length(qoss2[qoss2>suppred[qoss>cpptresh]])/length(qoss2)*100,"%",sep="",fill=T)
  cat("Percentage of points lying below the lower confidence limit=",length(qoss2[qoss2<infpred[qoss>cpptresh]])/length(qoss2)*100,"%",sep="",fill=T)
  }
if(NSeff==T)
  {
# If NSeff=T returns efficiency of the D-model alone
  if(is.null(qoss)==T)
  {
    cat("Error! Cannot compute efficiency if observed flow is not available. Exiting",fill=T)
    return()
  }
  qsimf$qtot[qsimf$qtot>(max(qoss)*2)]=max(qoss)*2
  eff1=1-sum((qsimf$qtot-qoss)^2)/sum((qoss-mean(qoss))^2)
  return(eff1)
  } else if(bluecat==T)
  {
  if(is.null(qoss)==F) 
    {
    eff=1-sum((medpred-qoss)^2)/sum((qoss-mean(qoss))^2)
    return(list(q_slow=qsimf$qtslow,q_quick=qsimf$qtquick,q_tot=qsimf$qtot,medpred=medpred,infpred=infpred,suppred=suppred,effsmodel=eff))
    } else
    {
    return(list(q_slow=qsimf$qtslow,q_quick=qsimf$qtquick,q_tot=qsimf$qtot,medpred=medpred,infpred=infpred,suppred=suppred))
    }
  }
  else
    {if(is.null(qoss)==F)
      {
      eff1=1-sum((qsimf$qtot-qoss)^2)/sum((qoss-mean(qoss))^2)
      return(list(q_slow=qsimf$qtslow,q_quick=qsimf$qtquick,q_tot=qsimf$qtot,effdmodel=eff1))
      } else
      {
      return(list(q_slow=qsimf$qtslow,q_quick=qsimf$qtquick,q_tot=qsimf$qtot))
      }
    }
}

#########################################################################
### Function to fit the PBF distribution on the stochastic simulated data
#########################################################################

fitPBF=function(paramd,ptot,kp,kptail)
{
#Computation or return period from k-moments by using default values for csi,zi,lambda
#Computation of lamdbda1 and lambdainf
lambda1k=(+1+(beta(1/paramd[2]/paramd[1]-1/paramd[2],1/paramd[2])/paramd[2])^paramd[2])^(1/paramd[2]/paramd[1])
lambda1t=1/(1-(+1+(beta(1/paramd[2]/paramd[1]-1/paramd[2],1/paramd[2])/paramd[2])^paramd[2])^(-1/paramd[2]/paramd[1]))
lambdainfk=gamma(1-paramd[1])^(1/paramd[1])
lambdainft=gamma(1+1/paramd[2])^(-paramd[2])
#Computation of return periods
Tfromkk=lambdainfk*ptot+lambda1k-lambdainfk
Tfromkt=lambdainft*ptot+lambda1t-lambdainft
Tfromdk=1/(1+paramd[1]*paramd[2]*((kp-paramd[4])/paramd[3])^paramd[2])^(-1/paramd[1]/paramd[2])
Tfromdt=1/(1-(1+paramd[1]*paramd[2]*((kptail-paramd[4])/paramd[3])^paramd[2])^(-1/paramd[1]/paramd[2]))
lsquares=sum(log(Tfromkk/Tfromdk)^2)+sum(log(Tfromkt/Tfromdt)^2)
if(is.na(lsquares)==T) lsquares=1E20
return(lsquares)
}
