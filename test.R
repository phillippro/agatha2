library(shiny)
library(magicaxis)
source('periodoframe.R')
source("periodograms.R")
source('functions.R',local=TRUE)
source('mcmc_func.R')
Nmax.plots <- 50
count0 <- 0
instruments <- c('HARPS','SOHPIE','HARPN','AAT','KECK','APF','PFS')
tol <- 1e-16
#trend <- FALSE
data.files <- list.files(path='data',full.name=FALSE)
tab <- read.table('data/ChallengeDataSet1_HARPS.dat',header=TRUE)
frange <- 1/10^c(4,0.1)
#rv.ls <- glst(t=tab[,1],y=tab[,2],err=tab[,3],ofac=1,fmin=frange[1],fmax=frange[2])
#Indices <- scale.proxy(tab[,4:ncol(tab),drop=FALSE])
Indices <- scale.proxy(tab[,6,drop=FALSE])
#Indices <- NULL
Nar <- 0
Nma <- 1
ofac <- 1
per.type <- per.type.seq <- 'MLP'
#basis <- 'linear'
basis <- 'natural'
SigType <- 'kepler'
#SigType <- 'circular'
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
renew <- FALSE
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

fout <- paste0('results/ChallengeSet1_',per.type,'_ofac',ofac,'_ARMA',Nma,Nar,'_NI',ncol(Indices),'_',SigType,'_',basis,'_fold',fold,'.pdf')
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))

if(per.type=='MLP') rv.ls <- MLP(t=tab[,1],y=tab[,2],dy=tab[,3],Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2],Indices=Indices,MLP.type='sub')
if(per.type=='BGLS') rv.ls <- bgls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='BFP'){
   rv.ls <- BFP(t=tab[,1],y=tab[,2],dy=tab[,3], Nma=Nma,Nar=Nar,model.type='man',Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify, renew=renew,progress=progress,noise.only=noise.only)
}
if(per.type=='GLS') rv.ls <- gls(t=tab[,1],y=tab[,2],err=tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='GLST') rv.ls <- glst(t=tab[,1],y=tab[,2],err=tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
if(per.type=='LS') rv.ls <- lsp(times=tab[,1],x=tab[,2],ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))

#if(SigType!='stochastic'){
tmp <- sigfit(per=rv.ls,data=tab,SigType=SigType,basis=basis,fold=fold)
rv.ls <- tmp$per
par.data <- addpar(c(),tmp$ParSig,1)

phase.plot(tmp,fold=fold)
  pp <- cbind(tmp$t,tmp$y,tmp$ey)
        colnames(pp) <- paste0(c('t','y','ey'),'_sig1')
        phase.data <- cbind(phase.data,pp)

        qq <- cbind(tmp$tsim,tmp$ysim)
        colnames(qq) <- paste0(c('tsim','ysim'),'_sig1')
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
Nsig.max <- 3
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
if(SigType!='stochastic'){
    source('additional_signals.R')
}
dev.off()
