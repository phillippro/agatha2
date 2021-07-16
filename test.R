library(shiny)
library(magicaxis)
source('periodoframe.R')
source("periodograms.R")
source('functions.R',local=TRUE)
Nmax.plots <- 50
count0 <- 0
instruments <- c('HARPS','SOHPIE','HARPN','AAT','KECK','APF','PFS')
tol <- 1e-16
#trend <- FALSE
data.files <- list.files(path='data',full.name=FALSE)
tab <- read.table('data/ChallengeDataSet2_HARPS.dat',header=TRUE)
out <- calcBF(data=tab,Nbasic=0,
            proxy.type='cum',
              Nma.max=1,Nar.max=1,
              Nproxy=3,progress=FALSE)


