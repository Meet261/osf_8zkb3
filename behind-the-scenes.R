##################################################################
####### Exploring sequential sampling models of behaviour ########
##################################################################

# Guy Hawkins, University of Amsterdam (guy.e.hawkins@gmail.com). 


### functions required for practical-main.R that user doesn't need to see ###
# load some packages
require(rtdists)
require(msm)
require(DEoptim)


# generate synthetic data from LBA or DDM
synthesise.data.1condition=function(parameters,n,model) {
  x=unlist(parameters) # for simple reference inside function
  if(length(n)==1) n=c(n,n)
  if(length(n)>2) stop("Incorrect cell size (n) specified. Specify length 1 or 2.")

  # object to hold synthetic data
  data=list(S=NULL,R=NULL,RT=NULL)
  nresps=2  # hard code for this application
  resps=1:nresps
  
  # simulate data from LBA
  if(model=="LBA") {
    x["s"]=1  # fix standard deviation of trial-to-trial variability of drift rates
    LBA.par.names=c("A","B","vc","ve","t0","s")
    if(!all(names(x)%in%LBA.par.names)) stop("Incorrect parameters specified for LBA.")
    # do some quick checks on parameter values
    if(x["A"]<0) stop("Maximum value of the start-point distribution must be positive (i.e., set A>0).")
    if(x["B"]<0) stop("Distance from the top of the start-point distribution to the threshold must be positive (i.e., set B>0).")
    if(x["vc"]<0 | x["ve"]<0) stop("Drift rates must be positive (i.e., set vc>0 and ve>0).")
    if(x["t0"]<0) stop("Non-decision time must be positive (i.e., set t0>0).")

    # if all parameters are ok, simulate data    
    b=x["A"]+x["B"]
    for(stim in resps) {
      v = rep(x["ve"],nresps)
      v[stim] = x["vc"] 
      # trials generated from LBA 
      tmp=rlba_norm(n=n[stim],b=b,A=x["A"],mean_v=v,sd_v=x["s"],t0=x["t0"])
      data$S=c(data$S,rep(stim,n[stim]))
      data$R=c(data$R,tmp$response)
      data$RT=c(data$RT,tmp$rt)
    }
    
  # simulate data from DDM
  } else if(model=="DDM") {
    x["d"]=0  
    DDM.par.names=c("a","z","v","t0","sv","sz","st0","d")
    if(!all(names(x)%in%DDM.par.names)) stop("Incorrect parameters specified for DDM.")
    # do some quick checks on parameter values
    if(x["a"]<0) stop("Response boundary must be positive (i.e., set a>0).")
    if(x["z"]<0 | x["z"]>1) stop("Relative start point (z) must be between 0 and 1.")
    if(x["t0"]<0) stop("Non-decision time must be positive (i.e., set t0>0).")
    if(x["sv"]<0) stop("SD for trial-to-trial variability in drift rates must be positive (i.e., set sv>0).")
    if(x["sz"]<0) stop("Range of uniform distribution on start point (sz) must be between 0 and 1.")
    if(x["st0"]<0) stop("Range of uniform distribution on non-decision time (st0) must be above 0 and does not generally exceed 0.3.")
    
    # if all parameters are ok, simulate data    
    for(stim in resps) {
      # upper boundary correct for stimulus 1 (positive v), lower boundary correct for stimulus 2 (negative v)
      v=x["v"]*c(1,-1)[stim]  
      # trials generated from DDM 
      tmp=rrd(n[stim],a=x["a"],v=v,t0=x["t0"],z=x["z"],d=x["d"],sz=x["sz"],sv=x["sv"],st0=x["st0"])
      # relevel so that upper boundary is level 1, rename levels
      tmp$response=relevel(tmp$response,ref="upper")
      levels(tmp$response)=c("1","2")
      data$S=c(data$S,rep(stim,n[stim]))
      data$R=c(data$R,as.numeric(as.character(tmp$response)))
      data$RT=c(data$RT,tmp$rt)
    }
  }
  
  # return synthetic data as data frame with stimulus (S), response (R), response time (RT)
  data=data.frame(data)
  data$S=factor(data$S)
  data$R=factor(data$R)
  # assign number of conditions to data frame
  attributes(data)$n.conditions=1
  data
}


# generate synthetic data from LBA or DDM for 2 conditions
synthesise.data.2condition=function(parameters,n,model) {
  n.pars=sapply(parameters,length)
  if(any(n.pars>2))   stop("Specify no more than two values for any parameter.")
  if(!any(n.pars==2)) stop("Specify at least one parameter to differ across the two conditions.")
  for(i in 1:length(parameters)) {
    if(length(parameters[[i]])==1) parameters[[i]]=rep(parameters[[i]],length.out=2)
  }
  # synthesise data for conditions 1 and 2
  tmp1=synthesise.data.1condition(parameters=sapply(sapply(parameters,"[",1),list),n=n,model=model)
  tmp2=synthesise.data.1condition(parameters=sapply(sapply(parameters,"[",2),list),n=n,model=model)
  # merge conditions 1 and 2 into a single data frame
  data=rbind(cbind(A=rep(1,nrow(tmp1)),tmp1), cbind(A=rep(2,nrow(tmp2)),tmp2))
  data$A=factor(data$A)
  # assign number of conditions to data frame
  attributes(data)$n.conditions=2
  data
}


# return summary statistics from data or model predictions
summary.statistics=function(data,print.quantiles=FALSE,quantiles=seq(.1,.9,.2)) {
  
  # return summary statistics for a single condition
  summarise.data.1condition=function(data,print.quantiles=FALSE,quantiles=seq(.1,.9,.2)) {
    corrects=data$S==data$R
    cat("\n Mean accuracy:",round(mean(corrects),3))
    cat("\n Mean correct response time:",round(mean(data$RT[corrects]),3),"s")
    cat("\n Mean error response time:",round(mean(data$RT[!corrects]),3),"s")
    if(print.quantiles) {
      cat("\n\n Correct response time quantiles:\n")
      print(round(quantile(data$RT[corrects],prob=quantiles),2))
      cat("\n Error response time quantiles:\n")
      print(round(quantile(data$RT[!corrects],prob=quantiles),2))
    }
  }

  if(attr(data,"n.conditions")==1) {
    summarise.data.1condition(data, print.quantiles=print.quantiles, quantiles=quantiles)
  } else {
    for(i in levels(data$A)) {
      cat(paste0("\n\nCondition A==",i,":\n"))
      summarise.data.1condition(subset(data,A==i), print.quantiles=print.quantiles, quantiles=quantiles)
    }
  }
}


# four different ways of plotting data
mean.plot=function(data) {
  corrects=data$S==data$R
  draw.ci=TRUE
  par(mar=c(5,5,2,1))
  # calculate mean and SE for correct and error RTs, and accuracy
  rt.mean=rev(tapply(data$RT,corrects,mean))
  rt.ci=rev(tapply(data$RT,corrects,function(x) qnorm(.975)*sd(x)/sqrt(length(x))))
  par(mar=c(5,5,2,1))
  ylims=range(c(rt.mean-rt.ci,rt.mean+rt.ci))
  plot(y=rt.mean,x=1:2,xaxt="n",ylim=ylims,xlim=c(1:2)+c(-1,1)*.2,
       main="Mean RT",xlab="",ylab="Response time (seconds)",pch=16,cex=1.5)
  axis(side=1,at=1:2,c("Correct","Error"))
  if(draw.ci) {
    arrows(x0=1:2,x1=1:2,y0=rt.mean-rt.ci,y1=rt.mean+rt.ci,
           code=3,angle=90,length=.1,lwd=1.5)
  }
}


histogram.plot=function(data) {
  # trim data in case there are any big RT outliers
  data=subset(data,RT<quantile(RT,probs=.99))
  corrects=data$S==data$R
  cor.hist=hist(data$RT[corrects],plot=F)
  # scale density of cor.hist by the proportion of correct responses
  cor.hist$density=cor.hist$density*mean(corrects)
  ymax=max(cor.hist$density)
  # only draw error histogram if there are enough errors in data
  plot.errors=FALSE
  if(sum(!corrects)>20) {
    plot.errors=TRUE
    err.hist=hist(data$RT[!corrects],plot=F)
    err.hist$density=err.hist$density*(1-mean(corrects))
    ymax=max(c(ymax,err.hist$density))
  } else {
    cat("\nInsufficient errors in data to generate a histogram of error response times.")
  }
  par(mar=c(5,5,2,1))
  # order the overlaying of corrects on errors and vice versa
  if(!plot.errors) {
    plot(cor.hist,freq=F,col=rgb(0,.8,0,.7),ylim=c(0,ymax),
         main="Histogram of RTs",xlab="Response time (seconds)")
  } else {
    if(mean(corrects)>.5) {  # overlay errors on corrects
      plot(cor.hist,freq=F,col=rgb(0,.8,0,.7),ylim=c(0,ymax),
           main="Histogram of RTs",xlab="Response time (seconds)")
      plot(err.hist,freq=F,col=rgb(1,0,0,.6),add=T)
    } else {  # overlay corrects on errors
      plot(err.hist,freq=F,col=rgb(1,0,0,.6),ylim=c(0,ymax),
           main="Histogram of RTs",xlab="Response time (seconds)")
      plot(cor.hist,freq=F,col=rgb(0,.8,0,.7),ylim=c(0,ymax),add=T)
    }
  }
}


get.quantiles.and.probs=function(data) {
  corrects=data$S==data$R
  p=mean(corrects)
  quantiles=seq(.1,.9,.2)
  # get quantiles for correct and error responses
  qs=cbind(quantile(data$RT[!corrects],probs=quantiles),quantile(data$RT[corrects],probs=quantiles))
  list(qs=qs,quantiles=quantiles,p=p)  
}


CDF.plot=function(data,model.predictions=NULL,plot.header="") {
  n.conditions=attr(data,"n.conditions")
  # Boolean to determine whether we overlay model predictions on data
  plot.model=ifelse(is.null(model.predictions),FALSE,TRUE)

  # plot data from single condition
  if(n.conditions==1) {
    if(plot.header=="") plot.header="Defective CDF"
    # compute quantiles
    tmp=get.quantiles.and.probs(data) ; data.qs=tmp$qs ; data.p=tmp$p ; quantiles=tmp$quantiles
    qs.range=range(data.qs)
    if(plot.model) { 
      tmp=get.quantiles.and.probs(model.predictions) ; model.qs=tmp$qs ; model.p=tmp$p 
      qs.range=range(c(qs.range,model.qs))       
    }
    xlims=qs.range+min(qs.range)*c(-1,1)*.05
    par(mar=c(5,5,2,1))
    matplot(x=data.qs,y=cbind(quantiles*(1-data.p),quantiles*data.p),type=ifelse(plot.model,"p","o"),
            pch=16,col=c("red2","green3"),cex=1.2,ylim=c(0,1),lty=1,lwd=2,xlim=xlims,
            xlab="Response time (seconds)",ylab="Probability of response",main=plot.header)
    # plot model predictions
    if(plot.model) {
      matpoints(x=model.qs,y=cbind(quantiles*c(1-model.p),quantiles*model.p),type="o",
             pch=4,col="black",cex=1.5,lty=1,lwd=2)
    }
  
  # plot data from multiple conditions
  } else if(n.conditions>1) {
    # compute quantiles
    data=by(data,data$A,get.quantiles.and.probs)
    qs.range=range(sapply(data,function(x) x$qs))
    if(plot.model) { 
      model=by(model.predictions,model.predictions$A,get.quantiles.and.probs)
      qs.range=range(c(qs.range,range(sapply(model,function(x) x$qs))))       
    }
    xlims=qs.range+min(qs.range)*c(-1,1)*.05
    par(mfrow=c(1,n.conditions),mar=c(4,4,2,1))
    # plot data
    for(i in 1:n.conditions) {
      matplot(x=data[[i]]$qs,y=cbind(data[[i]]$quantiles*(1-data[[i]]$p),data[[i]]$quantiles*data[[i]]$p),
              type=ifelse(plot.model,"p","o"),pch=16,col=c("red2","green3"),cex=1.2,ylim=c(0,1),lty=1,lwd=2,xlim=xlims,
              xlab="Response time (seconds)",ylab="Probability of response",main=paste(plot.header,"Condition",i))
      # plot model predictions
      if(plot.model) {
        matpoints(x=model[[i]]$qs,y=cbind(model[[i]]$quantiles*c(1-model[[i]]$p),model[[i]]$quantiles*model[[i]]$p),
                  type="o",pch=4,col="black",cex=1.5,lty=1,lwd=2)
      }
    }
  }
}


# quantile probability plot for 1 or N-condition experiment, with option to overlay model predictions
QP.plot=function(data,model.predictions=NULL,plot.header="") {
  n.conditions=attr(data,"n.conditions")
  if(plot.header=="") plot.header="Quantile probability plot"
  # Boolean to determine whether we overlay model predictions on data
  plot.model=ifelse(is.null(model.predictions),FALSE,TRUE)
  
  # plot data from single condition
  if(n.conditions==1) {
    # compute quantiles
    tmp=get.quantiles.and.probs(data) ; data.qs=tmp$qs ; data.p=tmp$p
    qs.range=range(data.qs)
    if(plot.model) { 
      tmp=get.quantiles.and.probs(model.predictions) ; model.qs=tmp$qs ; model.p=tmp$p 
      qs.range=range(c(qs.range,model.qs))       
    }
    ylims=qs.range+min(qs.range)*c(-1,1)*.05
    par(mar=c(5,5,2,1))
    # plot data
    plot(x=rep(c(1-data.p,data.p),each=nrow(data.qs)),y=data.qs,type="p",ylim=ylims,
         pch=16,col=rep(c("red2","green3"),each=nrow(data.qs)),cex=1.2,xlim=c(0,1),
         xlab="Probability of response",ylab="Response time (seconds)",main=plot.header)
    # plot model predictions
    if(plot.model) {
      for(i in 1:nrow(model.qs)) {
        points(x=c(1-model.p,model.p),y=model.qs[i,],type="o",pch=4,col="black",cex=1.5,lwd=1.5)
      }
    }
    
  # plot data from multiple conditions
  } else if(n.conditions>1) {
    # compute quantiles
    data=by(data,data$A,get.quantiles.and.probs)
    qs.range=range(sapply(data,function(x) x$qs))
    if(plot.model) { 
      model=by(model.predictions,model.predictions$A,get.quantiles.and.probs)
      qs.range=range(c(qs.range,range(sapply(model,function(x) x$qs))))       
    }
    ylims=qs.range+min(qs.range)*c(-1,1)*.05
    par(mar=c(5,5,2,1))
    # plot data
    plot(y=NULL,x=NULL,ylim=qs.range,xlim=c(0,1),
         xlab="Probability of response",ylab="Response time (seconds)",main=plot.header)
    for(i in 1:length(data)) {
      qs=data[[i]]$qs ; p=data[[i]]$p
      points(x=rep(c(1-p,p),each=nrow(qs)),y=qs,type="p",pch=16,col=rep(c("red2","green3"),each=nrow(qs)),cex=1.2)
    }
    # plot model predictions
    if(plot.model) {
      for(i in 1:nrow(model[[1]]$qs)) {
        ps=sapply(model,"[[","p")
        ps.ordered=c(rev(1-sort(ps)),sort(ps)) ; ord=names(ps.ordered)
        # order RT for each quantile. loop through error responses first, then corrects
        qs.ordered=NULL
        for(j in 1:length(model)) qs.ordered=c(qs.ordered,model[[ord[j]]]$qs[i,1])
        for(j in (length(model)+1):(2*length(model))) qs.ordered=c(qs.ordered,model[[ord[j]]]$qs[i,2])
        points(x=ps.ordered,y=qs.ordered,type="o",pch=4,col="black",cex=1.5,lwd=1.5)
      }
    }
  }
}


# plot data in various ways
plot.data=function(data,plot.type="mean",model.predictions=NULL,plot.header="") {
  n.conditions=attr(data,"n.conditions")
  # single condition
  if(n.conditions==1 & !(plot.type=="mean" | plot.type=="histogram" | plot.type=="CDF" | plot.type=="QP" | plot.type=="all")) {
      stop("Incorrect plot.type specified for one condition. Specify plot.type as \"mean\", \"histogram\", \"CDF\", \"QP\", or \"all\".") }

  # multiple conditions
  if(n.conditions>1 & !(plot.type=="CDF" | plot.type=="QP")) {
      stop("Incorrect plot.type specified for two conditions. Specify plot.type as \"CDF\", or \"QP\".") }

  # determine which plots to draw
  if(plot.type=="all") {
    # single condition - simultaneously display all four plot types
    par(mfrow=c(2,2))
    mean.plot(data)
    histogram.plot(data)
    CDF.plot(data,model.predictions=model.predictions,plot.header=plot.header)
    QP.plot(data,model.predictions=model.predictions,plot.header=plot.header)
    
  } else if(plot.type=="mean") {
    par(mfrow=c(1,1))
    mean.plot(data)
    
  } else if(plot.type=="histogram") {
    par(mfrow=c(1,1))
    histogram.plot(data)
      
  } else if(plot.type=="CDF") {
    par(mfrow=c(1,ifelse(n.conditions==1,1,2)))
    CDF.plot(data,model.predictions=model.predictions,plot.header=plot.header)
      
  } else if(plot.type=="QP") {
    par(mfrow=c(1,1))
    QP.plot(data,model.predictions=model.predictions,plot.header=plot.header)
  }  
}
  
  

### functions to perform parameter optimisation ###

# check whether drift rate or threshold is free to vary over conditions
check.fit.parameters=function(x) {
  if(sapply(x,length)["vc"]==2) { 
    parameterisation = "drift" 
  } else if(sapply(x,length)["B"]==2) { 
    parameterisation = "threshold" 
  } else { 
    stop("Specify two correct drift rates (vc) or two response thresholds (B)") 
  }
  parameterisation
}


# calculate log-likelihood of data given parameters. Assumes a single condition.
log.likelihood.1condition=function(x,data,model="LBA") {

  if(model=="LBA") {
    # get one parameter per accumulator 
    nresps = length(levels(data$R))
    resps = 1:nresps
    # Loop over stimuli and responses to get likelihoods
    out <- numeric(length(data$RT))
    for (stim in resps) {
      vs = rep(x["ve"],nresps)
      vs[stim] = x["vc"]
      for (resp in resps) {
        n1 = c(resp,resps[resps!=resp]) # index of 1st rate parameter
        isSR = data$S==levels(data$S)[stim] & data$R==levels(data$R)[resp]
        out[isSR] = n1PDF(t=data$RT[isSR],A=x["A"],b=x["b"],t0=as.numeric(x["t0"]),
                          mean_v=vs[n1],sd_v=x["sd_v"],silent=TRUE)
      }
    }
  }
  -sum(log(pmax(out,1e-10))) # For safety against log(0) problems
}


# calculate log-likelihood of data given parameters for a two-condition data set
# assume a drift rate OR threshold difference across conditions
get.log.likelihood=function(x,data,model="LBA",parameterisation="",par.names) {
  
  # call log.dens.like.1condition twice...
  # list to hold to-be-estimated parameters
  names(x)=par.names
  pars=list()  
  if(model=="LBA") {   # DDM not implemented
    fixedpar=1
    if(parameterisation=="one.condition") {
      pars=c(A=x["A"],b=x["A"]+x["B"],vc=x["vc"],ve=x["ve"],t0=x["t0"],sd_v=fixedpar)
      names(pars)=c("A","b","vc","ve","t0","sd_v")
    } else if(parameterisation=="drift") {
      pars[[1]]=c(A=x["A"],b=x["A"]+x["B"],vc=x["vc1"],ve=x["ve"],t0=x["t0"],sd_v=fixedpar)
      pars[[2]]=c(A=x["A"],b=x["A"]+x["B"],vc=x["vc2"],ve=x["ve"],t0=x["t0"],sd_v=fixedpar)
      names(pars[[1]])=names(pars[[2]])=c("A","b","vc","ve","t0","sd_v")
    } else if(parameterisation=="threshold") {
      pars[[1]]=c(A=x["A"],b=x["A"]+x["B1"],vc=x["vc"],ve=x["ve"],t0=x["t0"],sd_v=fixedpar)
      pars[[2]]=c(A=x["A"],b=x["A"]+x["B2"],vc=x["vc"],ve=x["ve"],t0=x["t0"],sd_v=fixedpar)
      names(pars[[1]])=names(pars[[2]])=c("A","b","vc","ve","t0","sd_v")
    }
  }
  
  # return negative log-likelihood for function minimisation
  if(parameterisation=="one.condition") {
    log.likelihood.1condition(x=pars,data=data,model=model)    
  } else if(parameterisation=="drift" | parameterisation=="threshold") {
    log.likelihood.1condition(x=pars[[1]],data=subset(data,A==1),model=model) + 
    log.likelihood.1condition(x=pars[[2]],data=subset(data,A==2),model=model)
  }
}


get.par.names=function(model="LBA",parameterisation="") {
  # get name of model parameters 
  list(LBA=list(one.condition=c("A","B","vc","ve","t0"),
                drift=c("A","B","vc1","vc2","ve","t0"),
            threshold=c("A","B1","B2","vc","ve","t0")),
       DDM=list(one.condition=c("a","z","v","t0","sv","sz","st0"),
                drift=c("a","z","v1","v2","t0","sv","sz","st0"),
            threshold=c("a1","a2","z","v","t0","sv","sz","st0")))[[model]][[parameterisation]]
}


fit.model=function(data,model="LBA",parameterisation="",max.iterations=250) {
  if(!(model=="LBA" | model=="DDM")) stop("Specify model as \"LBA\" r \"DDM\".")
  if(!(parameterisation=="one.condition" | parameterisation=="drift" | parameterisation=="threshold")) stop("Specify parameterisation as \"one.condition\" (i.e., estimate parameters from a single experimental condition), \"drift\" (i.e., estimate free drift rate parameter across two experimental conditions) or \n\"threshold\" (i.e., estimate free threshold parameter across two experimental conditions).")
  
  # get name of model parameters 
  par.names=get.par.names(model=model,parameterisation=parameterisation)
  # set lower and upper boundaries on optimisation routine  
  lowers=rep(0,length(par.names))
  uppers=rep(5,length(par.names))
  # estimate model parameters using differential evolution
  tmp=DEoptim(fn=get.log.likelihood,lower=lowers,upper=uppers,
              data=data,model=model,parameterisation=parameterisation,par.names=par.names,
              control=list(itermax=max.iterations))
  # assign model and parameterisation to fitted object
  attributes(tmp)$model=model
  attributes(tmp)$parameterisation=parameterisation
  tmp
}


synthesise.model.predictions=function(fitted.model) {
  model=attr(fitted.model,"model")
  parameterisation=attr(fitted.model,"parameterisation")
  n=10000   # number of trials per condition for simulated data
  x=fitted.model$optim$bestmem
  names(x)=NULL
  if(parameterisation=="one.condition") {
    pars=list(A=x[1],B=x[2],vc=x[3],ve=x[4],t0=x[5])
    synthesise.data.1condition(parameters=sapply(pars,list), n=n, model=model)
  } else if(parameterisation=="drift") {
    pars=list(A=x[1],B=x[2],vc=x[3:4],ve=x[5],t0=x[6])
    synthesise.data.2condition(parameters=pars, n=n, model=model)
  } else if(parameterisation=="threshold") {
    pars=list(A=x[1],B=x[2:3],vc=x[4],ve=x[5],t0=x[6])
    synthesise.data.2condition(parameters=pars, n=n, model=model)
  }
}


summarise.model.fit=function(fitted.model,data,print.fit.statistics=FALSE) {
  model=attr(fitted.model,"model")
  parameterisation=attr(fitted.model,"parameterisation")
  if(!(model=="LBA" | model=="DDM")) stop("Specify model as \"LBA\" r \"DDM\".")
  if(!(parameterisation=="one.condition" | parameterisation=="drift" | parameterisation=="threshold")) stop("Specify parameterisation as \"one.condition\" (i.e., estimate parameters from a single experimental condition), \"drift\" (i.e., estimate free drift rate parameter across two experimental conditions) or \n\"threshold\" (i.e., estimate free threshold parameter across two experimental conditions).")

  # get name of model parameters 
  pars = fitted.model$optim$bestmem
  names(pars) = get.par.names(model=model,parameterisation=parameterisation)

  cat("\n",model,"fit using",parameterisation,"parameterisation.\n")  
  cat("\n Parameter estimates:\n") ; print(round(pars,3))
  
  if(print.fit.statistics) {
    # maximised value of the log likelihood function
    ll=-fitted.model$optim$bestval
    # AIC: -2*(log likelihood) + 2*(number of free parameters)
    AIC=-2*ll+2*length(pars)
    # BIC: -2*(log likelihood) + (log of number of data points)*(number of free parameters)
    BIC=-2*ll+log(nrow(data))*length(pars)
    
    cat("\nFit statistics:")  
    cat("\n Log-likelihood:",round(ll,3))
    cat("\n AIC:",round(AIC,3))
    cat("\n BIC:",round(BIC,3))
  }
}
  
