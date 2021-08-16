library(shiny)
library(doMC)
library(foreach)
library(magicaxis)
source('periodoframe.R')
source("periodograms.R")
source('functions.R')
source('mcmc_func.R')
Nmax.plots <- 50
count0 <- 0
instruments <- c('HARPS','SOHPIE','HARPN','AAT','KECK','APF','PFS')
tol <- 1e-16
#trend <- FALSE
data.files <- list.files(path='data',full.name=FALSE)
#fin <- 'ChallengeDataSet1_HARPS.dat'
fin <- 'GJ667C_HARPS.dat'
tab <- read.table(paste0('data/',fin),header=TRUE)
frange <- 1/10^c(4,0.1)
var <- list()
var$files <- fin
#rv.ls <- glst(t=tab[,1],y=tab[,2],err=tab[,3],ofac=1,fmin=frange[1],fmax=frange[2])
#Indices <- scale.proxy(tab[,4:ncol(tab),drop=FALSE])
#Indices <- scale.proxy(tab[,6,drop=FALSE])
Indices <- NULL
Ncores <- 4
Niter <- 1e3
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
Nar <- 0
Nma <- 0
ofac <- 1
#mcf <- TRUE
mcf <- FALSE
per.type <- per.type.seq <- 'MLP'
#basis <- 'linear'
basis <- 'natural'
#SigType <- 'kepler'
SigType <- 'circular'
#SigType <- 'stochastic'
if(per.type!='BFP' & SigType=='stochastic'){
    SigType <- 'circular'
}
phase.data <- sim.data <- c()
#SigType <- 'circular'
if(SigType=='stochastic'){
    noise.only <- TRUE
    fold <- FALSE
}else{
    noise.only <- FALSE
    fold <- TRUE
}
if(per.type!='BFP' & per.type!='MLP'){
    Nma <- Nar <- 0
}
#renew <- FALSE
renew <- TRUE
progress <- FALSE
quantify <- FALSE
#noise.only <- FALSE

phase.plot <- function(df,fold=TRUE){
    plot(df$per$P,df$per$power,xlab='Period [day]',ylab='Power',log='x',type='l',main=per.type)
    abline(h=df$per$levels)
    if(fold){
        plot(df$t,df$y,xlab='Phase [day]',ylab='RV [m/s]',main=paste0(round(df$per$Popt,1),' d;RMS.white=',round(sd(tmp$res),1),';RMS.red=',round(sd(df$per$res),1)),col='white')
        lines(df$tsim,df$ysim,col='red')
    }else{
        plot(df$t,df$y,xlab='T-Tmin [day]',ylab='RV [m/s]',main=paste0(round(df$per$Popt,1),' d;RMS.white=',round(sd(tmp$res),1)),col='white')
#    lines(df$tsim,df$ysim,col=tcol('red',50))
        lines(df$tsim,df$ysim,col='red')
    }
    points(df$t,df$y)
#    plot(df$t,df$res,xlab='T',ylab='RV [m/s]',main=paste('RMS:',round(sd(df$y),1)))
}

fout <- paste0('results/ChallengeSet1_',per.type,'_ofac',ofac,'_ARMA',Nar,Nma,'_NI',ncol(Indices),'_',SigType,'_',basis,'_fold',fold,'_mc',mcf,'.pdf')
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))

if(per.type=='MLP') rv.ls <- MLP(t=tab[,1],y=tab[,2],dy=tab[,3],Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2],Indices=Indices,MLP.type='sub')
if(per.type=='BGLS') rv.ls <- bgls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='BFP'){
   rv.ls <- BFP(t=tab[,1],y=tab[,2],dy=tab[,3], Nma=Nma,Nar=Nar,model.type='man',Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify, renew=renew,progress=progress,noise.only=noise.only)
   if(FALSE){
   dev.new()
   plot(rv.ls$P,rv.ls$power,log='x',type='l')
   stop()
   }
}
if(per.type=='GLS') rv.ls <- gls(t=tab[,1],y=tab[,2],err=tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='GLST') rv.ls <- glst(t=tab[,1],y=tab[,2],err=tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='LS') rv.ls <- lsp(times=tab[,1],x=tab[,2],ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))

#if(SigType!='stochastic'){
tmp <- sigfit(per=rv.ls,data=tab,SigType=SigType,basis=basis,mcf=mcf,Niter=Niter)
#stop()
rv.ls <- tmp$per
ParSig <- par.data <- addpar(c(),tmp$ParSig,1)
#stop()
#if(SigType!='stochastic'){
    phase.plot(tmp,fold=TRUE)
#}else{
#    plot(tab[,1]-min(tab[,1]),tmp$yred+tmp$res,xlab='time',ylab='RV')
#    lines(tmp$tsim0,tmp$ysim.red,col='red')
#}
pp <- cbind(tmp$t,tmp$y,tmp$ysig0)
        colnames(pp) <- paste0(c('t','y','ysig'),'_sig1')
        phase.data <- cbind(phase.data,pp)

        qq <- cbind(tmp$tsim,tmp$ysim,tmp$ysim0)
        colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig1')
        sim.data <- cbind(sim.data,qq)
#}else{
#plot(df$per$P,df$per$power,xlab='Period [day]',ylab='Power',log='x',type='l',main=per.type)
#    abline(h=df$per$levels)
#}
cat('head(rv.ls$res)=',head(rv.ls$res),'\n')
cat('sd(rv.ls$res)=',sd(rv.ls$res),'\n')
#cat('Popt=',rv.ls$Popt,'\n')
#cat('sd(rv.ls$res)=',sd(rv.ls$res),'\n')
#cat('head(rv.ls$res)=',head(rv.ls$res),'\n')

##more signals
Nsig.max <- 2
if(SigType=='stochastic') Nsig.max <- 1
MLP.type <- 'sub'
per.data <- c()
Pmaxs <- c()
instrument <- 'HARPS'
ypar <- 'RV'
per.target <- 'DC1'
#Nma <- Nar <- 0
Inds <- 0
tits <- 'na'
fs <- 'lkjk'
sig.levels <- 3
cnames <- 'na'
name <- 'kjl'
ylabs <- 'y'
ylab <- 'y'
xlabs <- xlab <- 'x'
#source('additional_signals.R',local=TRUE)
t <- tab[,1]
dy <- tab[,3]
if(SigType!='stochastic' & Nsig.max>1){
    source('additional_signals.R')
    plot(rv.ls$P,rv.ls$power,log='x',type='l',main='signal 2')
}

res <- tmp$res
if(Nsig.max>1){
    ysig0 <- rowSums(phase.data[,paste0('ysig_sig',1:Nsig.max)])
    ysim <- rowSums(sim.data[,paste0('ysim0_sig',1:Nsig.max)])
}else{
    ysig0 <- phase.data[,'ysig_sig1']
    ysim <- sim.data[,'ysim0_sig1']
}
ysig <- ysig0+res
tsim0 <- tmp$tsim0
ParSig <- par.data

if(mcf & Nsig.max>1 & SigType!='stochastic'){
    fit <- mcfit(rv.ls,data=tab[,1:3],tsim=tsim0,Niter=Niter,SigType=SigType,basis=basis,ParSig=ParSig,Pconv=TRUE,Ncores=Ncores)
    ParSig <- fit$ParSig
    par.data <- fit$par.stat
    if(SigType=='stochastic'){
        ysig <- fit$yred
        ysim <- fit$ysim.red
    }else{
        ysig <- fit$ysig
        ysim <- fit$ysim.sig
    }
    res <- fit$res
}

plot(tab[,1],tab[,2],xlab='BJD',ylab='RV',main=paste('RMS =',round(sd(tab[,2]),1)))
plot(tab[,1]-min(tab[,1]),ysig,xlab='BJD',ylab='RV',main=paste('RMS =',round(sd(ysig),1)))
lines(tsim0,ysim,col='red')
plot(tab[,1]-min(tab[,1]),res,xlab='BJD',ylab='RV[m/s]',main=paste('RV residual; RMS =',round(sd(res),1)))

dev.off()
