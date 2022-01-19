# hymodbluecat
Software working in the R environment for applying the HyMod rainfall-runoff model to simulate/predict river flow data. The software provides two routines: one for parameter calibration (hymod.par) and one for producing the simulation/prediction with uncertainty estimation (hymod.sim). Confidence limits for the prediction are estimated by applying the BLUECAT method that transforms any deterministic model in a stochastic model, from which mean sotchastic simulation and confidence limits are obtained.
The following R libraries are needed to run the package:
- devtools
- DescTools
- DEoptim 

They can be installed in R with the following commands:

> install.packages("DescTools")

> install.packages("DEoptim")

> library(DescTools)

> library(DEoptim)

To install the software in R under the Linux operating system, please make sure that you have the libraries libxml2-dev and libcurl4-openssl-dev installed (additional libraries may be needed depending on the Linux version). Under Ubuntu 20.04 they can be installed by the following Linux command line:
> sudo apt-get install libxml2-dev libcurl4-openssl-dev

Then, "devtools" is needed into R, which can be installed by using the following command at the R console:

> install.packages("devtools")

> library(devtools)

Then, hymodbluecat is installed into R with the command:

> install_github("albertomontanari/hymodbluecat")

> library(hymodbluecat)

To install the software in R under the Windows operating system first download and install Rtools from http://cran.r-project.org/bin/windows/Rtools/) and then:

> install.packages("devtools")

> library(devtools)

> find_rtools()

> install_github("albertomontanari/hymodbluecat")

> library(hymodbluecat)

Please note that the latest version of R may be needed, so beware of the warnings that you may get during installation.

If you wish to reinstall the package, beware that you need to detach it first, with the instruction

> detach("package:hymodbluecat", unload=TRUE)

You may also need to restart R before reinstalling.

The software comes with two data sets described in Koutsoyiannis and Montanari (2021). They refer to the Arno and Sieve River Basins, in Italy.
To reproduce the case study of the Arno River basin as presented by Koutsoyiannis and Montanari (2021) the following R commands can be used:

> data(arnosubbiano)

> pr1=hymod.par(c(100,1,0.5,200,0.5),area=752,tdelta=86400,e=arnosubbiano[,3][1:7305],p=arnosubbiano[,2][1:7305],nstep=length(arnosubbiano[,2][1:7305]),qoss=arnosubbiano[,4][1:7305],qinitial=15,lower=c(10,0.1,0.1,0.1,0.1),upper=c(800,10,0.9,1000,100),opt="DEoptim",lambdaln=0.0001)

Before moving forward it is advisable to check that optimization gave back reasonable parameter values.

> pr2=hymod.sim(pr1$optim$bestmem,area=752,tdelta=86400,e=arnosubbiano[,3][1:7305],p=arnosubbiano[,2][1:7305],qinitial=15,qoss=arnosubbiano[,4][1:7305],resultcalib=pr1,bluecat=T,empquant=F,plot=T,m=100,predsmodel="mdn")

> pr3=hymod.sim(pr1$optim$bestmem,area=752,tdelta=86400,e=arnosubbiano[,3][7306:8036],p=arnosubbiano[,2][7306:8036],qinitial=15,qoss=arnosubbiano[,4][7306:8036],resultcalib=pr1,bluecat=T,empquant=F,plot=T,m=100,predsmodel="mdn")

To reproduce the case study of the Sieve River basin as presented by Koutsoyiannis and Montanari (2021) the following R commands can be used:

> data(sievefornacina)

> pr10=hymod.par(c(100,1,0.5,200,0.5),area=846,tdelta=3600,sievefornacina[3648:26304,3],sievefornacina[3648:26304,2],nstep=length(sievefornacina[3648:26304,3]),sievefornacina[3648:26304,4],qinitial=10,lower=c(10,0.0,0.1,0.1,0.1),upper=c(800,10,0.9,1000,100),opt="DEoptim",lambdaln=0)

> pr20=hymod.sim(pr10$optim$bestmem,area=846,tdelta=3600,e=sievefornacina[3648:26304,3],p=sievefornacina[3648:26304,2],qinitial=50,qoss=sievefornacina[3648:26304,4],resultcalib=pr10,bluecat=T,empquant=F,plot=T,predsmodel="mdn")

> pr30=hymod.sim(pr10$optim$bestmem,area=846,tdelta=3600,e=sievefornacina[(26305+3648):43848,3],p=sievefornacina[(26305+3648):43848,2],qinitial=50,qoss=sievefornacina[(26305+3648):43848,4],resultcalib=pr10,bluecat=T,empquant=F,plot=T,predsmodel="mdn")

A detailed explanation of the argument of the functions hymod.par and hymod.sim is given in the R help. To invoke it the following commands can be used:
> ?hymod.par

> ?hymod.sim

Please contact me if you would like additional help.



