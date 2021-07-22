getM0 <- function(e,omega,P,T,T0,type='primary'){
    Tp <- T-getphase(e,omega)*P
    ((T0-Tp)%%P)*2*pi/P
}
kep.mt2 <- function(m,e){
    tol = 1e-8
    E0 <- m
    Ntt <- 1e3
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
#        if(k==Ntt) cat('Keplerian solver does not converge:',e,m,E0,E1,'!\n')
        E0 <- E1
    }
    if(k==Ntt){
        cat('Keplerian solver does not converge!\n')
        cat('length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
    }
    return(E1)
}

nrc2 <- function(Nvar){
    nrow <- ceiling(sqrt(Nvar))
    if(Nvar<nrow*(nrow-1)){
        ncol <- nrow-1
    }else{
        ncol <- nrow
    }
    return(c(nrow,ncol))
}

nrc <- function(Nvar){
    nrow <- ceiling(Nvar/2)
    if(Nvar<2){
        ncol <- 1
    }else{
        ncol <- 2
    }
    return(c(nrow,ncol))
}

tv.per <- function(targets,ofac,data){
    for(target in targets){
        tab <- data[[target]]
        commandArgs <- function(trailingOnly=TRUE) c('NA',1000,100,ofac,'bgls','res')
        source('time_varying_periodogram.R',local=TRUE)
    }
}

addpar <- function(par.old,par.new,nsig){
    if(nsig==1){
        n0 <- names(par.new)
        if(any(n0=='A')|any(n0=='B')){
            names(par.new)[1:3] <- paste0(n0[1:3],1)
        }
        par <- par.new
    }else{
        n0 <- names(par.new)
        n1 <- n0
        par <- par.old
        if(any(n0=='A')|any(n0=='B')){
            n1[1:3] <- paste0(n0[1:3],nsig)
            names(par.new) <- n1
            Npar.noise <- length(par.new)-3
        }else if(any(grepl('omega',n0))){
            n1[1:5] <- gsub('1$',nsig,n0[1:5])
            names(par.new) <- n1
            Npar.noise <- length(par.new)-5
        }
        par <- c(par.old[-(length(par.old)-(Npar.noise:1)+1)],par.new)
    }
    return(par)
}

calc.1Dper <- function(Nmax.plots, vars,per.par,data){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    Nmas <- unlist(Nmas)
    Nars <- unlist(Nars)
    per.data <- c()
    phase.data <- c()
    sim.data <- c()
    par.data <- c()
    tits <- c()
    fs <- c()
    pars <- list()
    kk <- 1
    if(SigType=='stochastic' & per.type=='BFP'){
        noise.only <- TRUE
    }else{
        noise.only <- FALSE
    }
    MLP.type <- 'sub'
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds[[1]],Nma=Nmas[1],Nar=Nars[1])
                kk <- kk+1
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0,Nar=0)
                kk <- kk+1
            }
        }
    }
    Nvar <- min(length(pars),Nmax.plots)
    sig.levels <- c()
    cnames <- c()
    pers <- c()
    ylabs <- c()
    Pmaxs <- c()
    ypars <- c()
#    lapply(1:Nvar, function(i){
    for(i in 1:Nvar){
        var <- pars[[i]]$var
        if(length(per.target)==1){
            Nma <- as.integer(pars[[i]]$Nma)
            Nar <- as.integer(pars[[i]]$Nar)
            Inds <- as.integer(pars[[i]]$Inds)
        }
        per.type <- pars[[i]]$per.type
#        instrument <- paste(per.target,collapse='-')
        if(length(per.target)>1){
            instrument <- 'combined'
            subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
            tmp <- combine.data(data=subdata,Ninds=Inds,Nmas=Nmas,Nar=Nars)
            tab <- tmp$cdata
            idata <- tmp$idata
            colnames(tab) <- colnames(data[[1]])[1:3]
        }else{
            instrument <- per.target
            tab <- data[[per.target]]
        }
        ypar <- var
        ypars <- c(ypars,gsub(' ','',ypar))
        Indices <- NA
        if(ncol(tab)>3){
            Indices <- tab[,4:ncol(tab),drop=FALSE]
        }
        if(ypar==ns[1]){
            dy <- tab[,3]
        }else{
            dy <- rep(0.1,nrow(tab))
        }
        if(ypar!='Window Function'){
            y <- tab[,ypar]
            if(per.type=='GLST'){
                rv.ls <- glst(t=tab[,1],y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='GLS'){
                rv.ls <- gls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='BGLS'){
                rv.ls <- bgls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- expression('log(ML/'*ML[max]*')')
                name <- 'logML'
            }else if(per.type=='BFP'){
                if(exists('per.type.seq')){
                    if(per.type.seq=='BFP'){
                        quantify <- TRUE
                    }else{
                        quantify <- FALSE
                    }
                }else{
                    quantify <- FALSE
                }

###preselect Indices according to the value of Inds and Indices
                if(length(Inds)>0){
                    if(all(Inds==0)){
                        Indices <- NULL
                    }else{
                        Inds <- Inds[Inds>0]
                        Indices <- as.matrix(Indices[,Inds,drop=FALSE])
                        for(j in 1:ncol(Indices)){
                            Indices[,j] <- scale(Indices[,j])
                        }
                    }
                }else{
                    Indices <- NULL
                }
#                tmp <- c(Nma=Nma,Nar=Nar,model.type='man',Indices=NULL,
#                                                      ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify)
                rv.ls <- BFP(t=tab[,1],y=y,dy=dy, Nma=Nma,Nar=Nar,model.type='man',Indices=Indices, ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify, renew=renew,noise.only=noise.only)

###renew: every chi-square minimization start from the initial parameter values
                ylab <- 'ln(BF)'
                name <- 'logBF'
            }else if(per.type=='MLP'){
                if(length(per.target)>1){
                    rv.ls <- MLP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,Nma=0,Nar=0,Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],MLP.type=MLP.type)
                }else{
                    rv.ls <- MLP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],MLP.type=MLP.type)
                }
                ylab <- expression('log(ML/'*ML[max]*')')
                name <- 'logML'
            }else if(per.type=='LS'){
                rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=y,ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))
                ylab <- 'Power'
                name <- 'power'
            }
#            tit <- paste('Periodogram:',per.type,'; Target:',instrument,'; Observable',ypar)
            tit <- paste0(per.type,'; ',instrument,';', ypar,';1 signal')
            if(!exists('Nma')){
                Nma <- 0
            }
            if(!exists('Nna')){
                Nna <- 0
            }
            if(!exists('Inds')){
                Inds <- 0
            }
#            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',per.type,'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',per.type,'_AR',paste(Nars,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
        }else{
            rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=rep(1,nrow(tab)),ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))
            tit <- paste0('LS;',instrument,';',ypar)
            pt <- 'LS'
            if(!exists('Nma')){
                Nma <- 0
            }

            if(!exists('Nar')){
                Nar <- 0
            }
            if(!exists('Inds')){
                Inds <- 0
            }
#           f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',pt,'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',pt,'_AR',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            ylab <- 'Power'
            name <- 'power'
        }
        ylabs <- c(ylabs,ylab)
        tits <- c(tits,tit)
        fs <- c(fs,f)
        pers <- c(pers,per.type)
###plot
#        plotname <- paste("plot", i, sep="")
        if(per.type=='MLP' | per.type=='BGLS'){
            yy  <- rv.ls$power-max(rv.ls$power)
            rv.ls$sig.level <- NULL#max(yy)-log(c(10,100,1000))
        }else{
            yy <- rv.ls$power
        }
        if(!is.null(per.data)){
            if(nrow(per.data)>length(yy)){
                rv.ls$P <- c(rv.ls$P,rv.ls$P[length(rv.ls$P)])
                yy <- c(yy,yy[length(yy)])
            }else if(nrow(per.data)<length(yy)){
                rv.ls$P <- rv.ls$P[-length(rv.ls$P)]
                yy <- yy[-length(yy)]
            }
        }
        if(length(rv.ls$sig.level)<3){
            sig.levels <- cbind(sig.levels,c(rv.ls$sig.level,rep(NA,3-length(rv.ls$sig.level))))
        }else{
            sig.levels <- cbind(sig.levels,rv.ls$sig.level)
        }
        if(i==1) per.data <- cbind(per.data,rv.ls$P)
        per.data <- cbind(per.data,yy)
        Pmaxs <- c(Pmaxs,format(per.data[which.max(yy),1],digit=2))
        inds <- (ncol(per.data)-1):ncol(per.data)
        if(i==1)  cnames <- c(cnames,'P')
        cnames <- c(cnames,paste0(pers[i],'1signal:',gsub(' .+','',ypar),':',name))

####calculate the Keplerian fit
        tmp <- sigfit(per=rv.ls,data=data,SigType=SigType)
###update the output from periodogram
        rv.ls <- tmp$per

        pp <- cbind(tmp$t,tmp$y,tmp$ysig0)
        colnames(pp) <- paste0(c('t','y','ysig'),'_sig1')
        phase.data <- cbind(phase.data,pp)

        qq <- cbind(tmp$tsim,tmp$ysim,tmp$ysim0)
        colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig1')
        sim.data <- cbind(sim.data,qq)
        par.data <- par.data <- addpar(c(),tmp$ParSig,1)

    if(Nsig.max>1){
            if(per.type==per.type.seq){
                if(length(per.target)>1){
                    Nma <- 0
                    Nar <- 0
                    Inds <- 0
                }
                source('additional_signals.R',local=TRUE)
            }
        }
###attach common data
        phase.attach <- cbind(tab[,1],y,dy,tmp$res)
        colnames(phase.attach) <- c('t0','y0','ey0','res')
        sim.attach <- t(t(tmp$tsim0))
        colnames(sim.attach) <- 'tsim0'
        phase.data <- cbind(phase.data,phase.attach)
        sim.data <- cbind(sim.data,sim.attach)
    }
    colnames(per.data) <- cnames
    if(!exists('Nsig.max')){
        Nsig.max <- 1
    }
    if(!exists('Nma')){
        Nma <- 0
    }
    if(!exists('Inds')){
        Inds <- 0
    }
    fname <- paste0(paste(per.target,collapse='_'),'_',paste(ypars,collapse='.'),'_',paste(per.type,collapse=''),'_MA',paste(Nma,collapse=''),'proxy',paste(Inds,collapse='.'),'_',Nsig.max,'sig_',paste(Pmaxs,collapse='d'),'d')
    fname <- paste0(paste(per.target,collapse='_'),'_',paste(ypars,collapse='.'),'_',paste(per.type,collapse=''),'_AR',paste(Nma,collapse=''),'proxy',paste(Inds,collapse='.'),'_',Nsig.max,'sig_',paste(Pmaxs,collapse='d'),'d')
    return(list(per.data=per.data,phase.data=phase.data,sim.data=sim.data,par.data=par.data,tits=tits,pers=pers,levels=sig.levels,ylabs=ylabs,fname=fname,fs=fs))
}

sigfit <- function(per,data,SigType='circular',basis='natural',fold=TRUE){
###This function is to modify the output of various periodograms to give residual, model prediction, and optimal parameters as well as posterior/likelihood samples
##x is a list
##SigType is either circular or kepler
    ParSig <- par.opt <- unlist(per$par.opt)
    par.list <- as.list(per$par.opt)

    if(any(names(per)=='data')){
        data <- per$data
    }else{
        per$data <- data
    }
    if(!any(names(per$df)=='data')){
        per$df$data <- data
    }
    tmin <- min(data[,1])%%2400000
    t <- data[,1]%%2400000-tmin
    tsim <- seq(0,max(t),length.out=1e4)
    if(!any(names(per)=='ysims')) per$ysims <- ysims <- 0
    popt <- per$Popt
    save.data <- FALSE

    if(SigType=='circular'){
        ysim.sig <- par.opt['A']*cos(2*pi/popt*tsim)+par.opt['B']*sin(2*pi/popt*tsim)
        ysim.all <- ysim.sig
        if(any(names(par.opt)=='gamma')) ysim.all <- ysim.all+par.opt['gamma']
        if(any(names(par.opt)=='beta')) ysim.all <- ysim.all+par.opt['beta']*tsim
        ysig0 <- par.opt['A']*cos(2*pi/popt*t)+par.opt['B']*sin(2*pi/popt*t)
        ysig <- per$res+ysig0
        res <- per$res
        if(any(names(per)=='df')){
            df <- per$df
            if(any(names(df)=='NI') & any(names(df)=='data') & any(names(df)=='Indices')){
                if(!any(names(df)=='Nma')) df$Nma <- 0
                if(!any(names(df)=='NI')) df$NI <- 0
                if(!any(names(df)=='Nar')) df$Nar <- 0
                sim <- CircularSim(par.list,df,popt,tsim)
                ysim.sig <- sim$ysig
                ysim.all <- sim$y
                fit <- CircularSig(par.list,df)
                ysig <- fit$res+fit$ysig
                ysig0 <- fit$ysig
                res <- fit$res
            }
        }
        ParSig <- c(P=popt,ParSig)
    }else if(SigType=='kepler'){
###Keplerian fitting
        df <- per$df
        if(!any(names(df)=='Nma')) df$Nma <- 0
        if(!any(names(df)=='NI')) df$NI <- 0
        if(!any(names(df)=='Nar')) df$Nar <- 0
        fit <- KeplerFit(per,basis=basis)
        sig <- KeplerSig(fit$ParKep,df,basis=basis)
        sim <- KeplerSim(fit$ParKep,df,tsim,basis=basis)#with trend
        res <- sig$res
        ParSig <- unlist(fit$ParKep)#unlist
        ysig0 <- sig$ysig
        ysig <- sig$res+ysig0
        per$res <- per$data[,2]-sig$ysig
        popt <- fit$ParKep$P1
        ysim.all <- sim$y
        ysim.sig <- sim$ysig
        per$ysims <- per$ysims+ysim.sig
    }else if(SigType=='stochastic'){
###Stochastic fitting
        df <- per$df
        fit <- CircularSig(par.list,df)
        ysig0 <- fit$yred
        ysig <- fit$res+ysig0
        res <- per$res <- fit$res
        sim <- CircularSim(par.list,df,popt,tsim)
        ysim.sig <- sim$yred
        ysim.all <- sim$y
        per$ysims <- per$ysims+ysim.sig
    }
####output
    tsim1 <- tsim%%popt
    inds <- sort(tsim1,index.return=TRUE)$ix
    tsim2 <- tsim1[inds]
    ysim2 <- ysim.sig[inds]
    tsim0 <- tsim
    ysim0 <- ysim.sig
    if(fold){
        ts <- t%%popt
        tsims <- tsim2
        ysims <- ysim2
    }else{
        ts <- t
        tsims <- tsim
        ysims <- ysim.sig
    }
    return(list(per=per,t=ts,y=as.numeric(ysig),ey=data[,3],res=as.numeric(res),ysig0=ysig0,tsim0=tsim0,ysim0=ysim0,tsim=tsims,ysim=ysims,ParSig=ParSig))
}

phase1D.plot <- function(phase.data,sim.data,tits,download=FALSE,index=NULL,repar=TRUE){
    if(repar){
        if(is.null(index)){
            par(mfrow=c(ceiling(Nmax.plots/2),2),mar=c(5,5,3,1),cex.lab=1.5,cex.axis=1.5,cex=1,cex.main=1.0)
        }
        if(download & is.null(index)){
            par(mfrow=c(2,2),mar=c(5,5,3,1),cex.lab=1.2,cex.main=0.8,cex.axis=1.2,cex=1)
        }
    }
    if(!is.null(index)){
        inds <- index
        titles <- rep('',inds)
    }else{
        inds <- 1:length(grep('y_sig',colnames(phase.data)))
        titles <- tits
    }
    ysig0 <- 0
    ysim0 <- 0
    for(i in inds){
        t <- phase.data[,paste0('t_sig',i)]
        y <- phase.data[,paste0('y_sig',i)]
        ey <- phase.data[,'ey0']

        tsim <- sim.data[,paste0('tsim_sig',i)]
        ysim <- sim.data[,paste0('ysim_sig',i)]

        ysig0 <- ysig0+y
        ysim0 <- ysim0+sim.data[,paste0('ysim0_sig',i)]

        ylab <- 'y'
        plot(t,y,xlab='Phase [day]',ylab='Y',main=titles[i],col='white')
#    lines(df$tsim,df$ysim,col=tcol('red',50))
        lines(tsim,ysim,col='red',lwd=3)
        points(t,y,col='black')
        arrows(t,y-ey,t,y+ey,length=0.03,angle=90,code=3,col=tcol('black',50))
    }
###total fit
    t <- phase.data[,'t0']
    ey <- phase.data[,'ey0']
    res <- phase.data[,'res']
    y <- res+ysig0
    tsim0 <- sim.data[,'tsim0']

    plot(t,y,xlab='BJD [day]',ylab='Y',main=gsub('\\d signal','combined fit',titles[i]),col='white')
    lines(tsim0+min(t),ysim0,col='red',lwd=3)
    points(t,y,col='black')
    try(arrows(t,y-ey,t,y+ey,length=0.03,angle=90,code=3,col=tcol('black',50)),TRUE)
    legend('topright',bty='n',legend=paste('RMS =',round(sd(y),1),'\n'))

###residual
    plot(t,res,xlab='BJD [day]',ylab='Y',main=gsub('\\d signal','residual',titles[i]))
    try(arrows(t,res-ey,t,res+ey,length=0.03,angle=90,code=3,col=tcol('black',50)),TRUE)
    legend('topright',bty='n',legend=paste('RMS =',round(sd(res),1),'\n'))
}

combined.plot <- function(per.data,phase.data,sim.data,tits,pers,levels,ylabs,SigType='circular',download=FALSE,index=NULL){
    per1D.plot(per.data,tits,pers,levels,ylabs,download=download,index=index,SigType=SigType)
    phase1D.plot(phase.data,sim.data,tits=tits,download=download,index=index,repar=FALSE)
}

per1D.plot <- function(per.data,tits,pers,levels,ylabs,download=FALSE,index=NULL,SigType='circular'){
    if(is.null(index)){
        par(mfrow=c(ceiling(Nmax.plots/2),2),mar=c(5,5,3,1),cex.lab=1.5,cex.axis=1.5,cex=1,cex.main=1.0)
    }
    if(download & is.null(index)){
        par(mfrow=c(2,2),mar=c(5,5,3,1),cex.lab=1.2,cex.main=0.8,cex.axis=1.2,cex=1)
    }
    P <- per.data[,1]
    if(!is.null(index)){
        inds <- index
        titles <- rep('',inds)
    }else{
        inds <- 1:(ncol(per.data)-1)
        titles <- tits
    }
    for(i in inds){
        power <- per.data[,i+1]
        ylab <- ylabs[i]#gsub('.+:','',colnames(per.data)[i+1])
        per.type <- gsub('[[:digit:]]signal:.+','',colnames(per.data)[i+1])
        f1 <- gsub('signal:.+','',colnames(per.data)[i+1])
        Nsig <- gsub('[A-Z]','',f1)
        if(pers=='BGLS'|pers=='MLP'|(pers=='BFP' & SigType!='stochastic')){
            ymin <- median(power)
        }else{
            ymin <- max(0,min(power))
        }
        if(grepl('Window',titles[i])){
            ylim <- c(ymin,max(power)+0.15*(max(power)-ymin))
        }else{
#            ylim <- c(ymin,max(max(power)+0.15*(max(power)-ymin),levels[which(!is.na(levels[,i])),i]))
            ylim <- c(ymin,max(power)+0.15*(max(power)-ymin))
        }
        plot(P,power,xlab='Period [day]',ylab=ylab,xaxt='n',log='x',type='l',main=titles[i], ylim=ylim)
        magaxis(side=1,tcl=-0.5)
        if(!grepl('Window',titles[i])){
            abline(h=levels[,i],lty=2)
        }
        p <- show.peaks(ps=P,powers=power,levels=levels[,i])
        if(!is.matrix(p)){
            pmaxs <- p[1]
            power.max <- p[2]
        }else{
            pmaxs <- p[,1]
            power.max <- p[,2]
            if(length(pmaxs)>4){
                pmaxs <- p[1:2,1]
                power.max <- p[1:2,2]
            }
        }
        if(length(pmaxs)>0){
            par(xpd=TRUE)
            offset <- c(0.08*(max(power)-ymin), 0.02*(max(power)-ymin))
            pmaxs <- pmaxs[1]
            power.max <- power.max[1]
            text(pmaxs,power.max+offset[1],pos=3,labels=format(pmaxs,digit=4),col='red',cex=1.0)
            try(arrows(pmaxs,power.max+offset[1],pmaxs,power.max+offset[2],col='red',length=0.05),TRUE)
            par(xpd=FALSE)
        }
    }
}

per2D.data <- function(vars,per.par,data){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    Nmas <- unlist(Nmas)
    pars <- list()
    kk <- 1
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds[[1]],Nma=Nmas[1])
                kk <- kk+1
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0)
                kk <- kk+1
            }
        }
    }
    i <- 1
    if(length(per.target)==1){
        Nma <- as.integer(pars[[i]]$Nma)
        Inds <- as.integer(pars[[i]]$Inds)
    }
    Indices <- NA
    per.type <- pars[[i]]$per.type
    var <- pars[[i]]$var
    if(length(per.target)>1){
        instrument <- 'combined'
        subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
        tmp <- combine.data(data=subdata,Ninds=Inds,Nmas=Nmas)
        tab <- tmp$cdata
        idata <- tmp$idata
        colnames(tab) <- colnames(data[[1]])[1:3]
    }else{
        instrument <- per.target
        tab <- data[[per.target]]
        if(ncol(tab)>3){
            Indices <- as.matrix(tab[,4:ncol(tab),drop=FALSE])
        }
    }
    ypar <- var
    t <- tab[,1]%%2400000#min(tab[,1])
    y <- tab[,yvar]
    dy <- tab[,3]
    if(length(per.target)==1){
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=Nma,Inds=Inds,Indices=Indices)
    }else{
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=0,Inds=0,Indices=Indices)
    }
    x2 <- mp$tmid
    y2 <- mp$P
    z2 <- mp$powers
    z2.rel <- mp$rel.powers
    fname <- paste0(paste(per.target,collapse='_'),'_MP_',paste(per.type,collapse=''),'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'))
    if(length(per.target)==1){
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,fname=fname,ypar=ypar))
    }else{
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,subdata=subdata,idata=idata,fname=fname,ypar=ypar))
    }
}

plotMP <- function(vals,pars){
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
    if(length(per.target)>1){
        subdata <- vals$subdata
        idata <- vals$idata
    }
    ypar <- vals$ypar
    t <- vals$t
    y <- vals$y
    dy <- vals$dy
    xx <- vals$xx
    yy <- vals$yy
    zz <- vals$zz
    zz.rel <- vals$zz.rel
    source('MP_plot.R',local=TRUE)
}

calcBF <- function(data,Nbasic,proxy.type,Nma.max,Nar.max,groups=NULL,Nproxy=NULL,Npoly=c(2,0),progress=FALSE){
##add Nar.max
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    NI.max <- ncol(data)-3
    NI.inds <- list(0)
    if(NI.max>0){
        NI.inds <- list()
        Nvary <- NI.max-Nbasic
        if(proxy.type=='cum' & Nvary>0 & Nproxy>0){
            NI.inds[[1]] <- Nbasic
            Indices <- data[,4:ncol(data),drop=FALSE]
            cors <- c()
            ###detrend the data first
            if(Npoly[1]>0){
                x <- t
                p <- lm(y~poly(x,Npoly[1]))
                y1 <- residuals(p)
            }else{
                y1 <- y
            }
            for(j in 1:ncol(Indices)){
                if(sd(Indices[,j])==0){
                    cors <- c(cors,0)
                }else{
                    if(Npoly[2]>0){
                        z <- Indices[,j]
                        p <- lm(z~poly(x,Npoly[2]))
                        z1 <- residuals(p)
                    }else{
                        z1 <- Indices[,j]
                    }
                    cors <- c(cors,abs(cor(z1,y1)))
                }
            }
            inds <- sort(cors,decreasing=TRUE,index.return=TRUE)$ix
            if(Nproxy>Nbasic){
                for(j in 1:(Nproxy-Nbasic)){
                    NI.inds[[j+1]] <- inds[1:j]
                }
            }
        }else if(proxy.type=='group' & Nvary>0){
            NI.inds <- lapply(1:(length(groups)+1),function(i) NI.inds[[i]] <- list())
            if(Nbasic>0){
                NI.inds[[1]] <- 1:Nbasic
            }else{
                NI.inds[[1]] <- 0
            }
            groups <- sort(as.integer(groups))
            for(j in 1:length(groups)){
                if(j==1){
                    NI.inds[[j+1]] <- 1:groups[j]
                }else{
                    if(Nbasic>0){
                        NI.inds[[j+1]] <- c(1:Nbasic,(groups[j-1]+1):groups[j])
                    }else{
                        NI.inds[[j+1]] <- (groups[j-1]+1):groups[j]
                    }
                }
            }
        }else if(proxy.type=='man'){
            NI.inds <- groups
        }else{
            NI.inds <- list(list(Nbasic:NI.max))
        }
    }
    Nmas <- 0:Nma.max
    Nars <- 0:Nar.max
    if(ncol(data)>3){
#        out <- BFP.comp(data, Nmas=0:Nma.max,Nars=0:Nar.max,NI.inds=NI.inds,progress=progress)
        out <- bfp.inf.progress(data,Nmas=Nmas,Nars=Nars,NI.inds=NI.inds)
    }else{
        out <- bfp.inf.progress(data,Nmas=Nmas,Nars=Nars,NI.inds=0)
#        out <- BFP.comp(data, Nmas=0:Nma.max,Nars=0:Nar.max,NI.inds=0,progress=progress)
    }
#    out$logBF
#    if(!is.matrix(out$logBFs)){
    lnBFs <- flatten2d(out$lnBF)
#    }
    return(list(Inds=NI.inds,Inds.opt=out$NI.opt,Nars=0:Nar.max,Nmas=0:Nma.max,Nma.opt=out$Nma.opt,Nar.opt=out$Nar.opt,lnBF=lnBFs,extra=out))
#    return(out)
}

flatten3d <- function(arr){
    dn <- dimnames(arr)
    coln <- gsub('\\d','',c(dn[[1]][1],dn[[2]][1],dn[[3]][1]))
    nn <- outer(outer(gsub('[a-z]|[A-Z]','',dn[[1]]),gsub('[a-z]|[A-Z]','',dn[[2]]),paste),gsub('[a-z]|[A-Z]','',dn[[3]]),paste)
    ns <- t(sapply(1:length(nn),function(i) unlist(strsplit(nn[i],' '))))
    tmp <- data.frame(cbind(ns,flatten(arr)))
    colnames(tmp) <- c(coln,'val')
#    tmp[,1:3] <- gsub('[a-z]|[A-Z]','',tmp[,1:3])
    tmp
}

flatten2d <- function(arr){
    Ncol <- dim(arr)[2]*dim(arr)[3]
    Nrow <- dim(arr)[1]
    dn <- dimnames(arr)
    cn  <- paste0('ARMA(',gsub(' ',',',outer(gsub('[a-z]|[A-Z]','',dn[[3]]),gsub('[a-z]|[A-Z]','',dn[[2]]),paste)),')')
    out <- array(NA,dim=c(Nrow,Ncol))
    colnames(out) <- cn
    rownames(out) <- unlist(dimnames(arr)[1])
    j <- 1
    for(i in 1:Nrow){
        out[i,] <- flatten(arr[i,,])
    }
    out
}

MCMC.panel <- function(){
    id <- 'HD020794_TERRA_1AP1_ervab6ap_ccf'
    Niter <- 1.0e3
    Nbin.per <- 1
    nbin.per <- 1
    tem <- 1
    inicov <- 1e-3
    Pini <- 200#day
    noise.model <- 'ARMA05'#noise.model: white, GP(R), ARMA, TJ(Ntj=1,noise vary with RHK or SA index, the third column of HARPS data), TJ(Ntj=3,vary with FWHM, BIS, RHK), TARMA, TGP, ARMATJ(ARMA+TJ), GPTJ,TJAR(model the RV contributed by index as a AR(p)-like model), ARMATJAR(ARMA+TJAR), PSID (previous-subsequent index dependent, this model is similar to TJAR but without time-varying/index-dependent jitter), ARMAPSID(ARMA+PSID)
    period.par <- 'logp'
    Ncores <- 1
    Np <- 1
    mode <- 'data'#data,sim
    Dtye <- 'D'#Dtype: DE:differential exclusing the target aperture, D: different including all aperture, N: no dependence on differential RV
    Nw <- 1#fit models to multiple wavelength data sets simultaneously
    prior.type <- 'mt'
    calibration <- 0#
    commandArgs <- function(trailingOnly=TRUE){
        cat('args=',c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
           ,'\n')
        c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
    }
    source('../mcmc_red.R',local=TRUE)
    return(list(folder=folder,pdf=gsub('.+/','',pdf.name)))
}
data.distr <- function(x,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    xs = sort(x)
    x1per = max(min(xs),xs[floor(length(xs)*0.01)])
    x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    return(c(x1per=x1per,x99per=x99per,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
}
show.peaks <- function(ps,powers,levels=NULL,Nmax=5){
    if(is.null(levels)) levels <- max(max(powers)-log(150),median(powers))
    ind <- which(powers==max(powers) | (powers>(max(powers)-log(100)) & powers>max(levels)))
    if(max(powers)-min(powers)<5) ind <- which.max(powers)
    pmax <- ps[ind]
    ppmax <- powers[ind]
    j0 <- 1
    p0 <- pmax[1]
    pp0 <- ppmax[1]
    pms <- p0
    pos <- pp0
    if(length(pmax)>1){
        for(j in 2:length(pmax)){
            if(abs(pmax[j]-p0) < 0.1*p0){
                if(ppmax[j]>pp0){
                    j0 <- j
                    p0 <- pmax[j]
                    pp0 <- ppmax[j0]
                    pms[length(pms)] <- p0
                    pos[length(pos)] <- pp0
                }
			    }else{
                j0 <- j
                p0 <- pmax[j]
                pp0 <- ppmax[j0]
                pms <- c(pms,p0)
                pos <- c(pos,pp0)
            }
        }
    }else{
        pms <- pmax
        pos <- ppmax
    }
    if(length(pms)>Nmax){
      pms <- pms[1:Nmax]
      pos <- pos[1:Nmax]
    }
    return(cbind(pms,pos))
}
