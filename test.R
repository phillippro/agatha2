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
rv.ls <- glst(t=tab[,1],y=tab[,2],err=tab[,3],ofac=1,fmin=frange[1],fmax=frange[2])
#SigType <- 'kepler'
SigType <- 'circular'

if(SigType=='kepler'){
tmp <- sigfit(per=rv.ls,SigType=SigType)
rv.ls <- tmp$per
}
cat('head(rv.ls$res)=',head(rv.ls$res),'\n')
cat('sd(rv.ls$res)=',sd(rv.ls$res),'\n')
cat('Popt=',rv.ls$Popt,'\n')
#cat('sd(rv.ls$res)=',sd(rv.ls$res),'\n')
#cat('head(rv.ls$res)=',head(rv.ls$res),'\n')

##more signals
Nsig.max <- 3
ofac <- 1
per.type <- per.type.seq <- 'GLST'
per.data <- c()
Pmaxs <- c()
instrument <- 'HARPS'
ypar <- 'RV'
per.target <- 'DC1'
Nma <- Nar <- 0
Inds <- 0
tits <- 'na'
fs <- 'lkjk'
sig.levels <- 3
cnames <- 'na'
name <- 'kjl'
ylabs <- 'y'
ylab <- 'y'
xlabs <- xlab <- 'x'
source('additional_signals.R',local=TRUE)

# rv.ls <- BFP(t=tab[,1],y=tab[,2],dy=tab[,3],                           Nma=0,Nar=0,model.type='man',Indices=NA,                            ofac=0.1,fmin=,fmax=frange[2],quantify=quantify, renew=renew)

