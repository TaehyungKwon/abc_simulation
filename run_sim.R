rm(list=ls())
####################################################################################
############################MADE BY TAEHYUNG KWON Nov 11 2019#######################
####################################################################################
genotype=function(geno){
	geno.vec=as.vector(geno)
	tn.size=length(geno.vec)
	ext.allele.vec=c(ceiling(geno.vec/2), floor(geno.vec/2))
	t.idx.allele=sample(x=c(0:1), size=tn.size, replace=T)
	allele.vec=ext.allele.vec[(1:tn.size) + tn.size * t.idx.allele]
	return(allele.vec)
}

evolve_first=function(pop,FM0,PZ,MNS,ZMS){
	pop.m=pop[which(pop[,4]==1),]
	pop.f=pop[which(pop[,4]==0),]
	fit.m=pop.m[,5]/sum(pop.m[,5])
	fit.f=pop.f[,5]/sum(pop.f[,5])
	
	parents=rbind(cbind(pop.m[sample(1:nrow(pop.m),size=nM,prob=fit.m,replace=T),],pop.f[sample(1:nrow(pop.f),size=nM,prob=fit.f,replace=T),]),cbind(pop.m[sample(nrow(pop.m),size=nF,prob=fit.m,replace=T),],pop.f[sample(nrow(pop.f),size=nF,prob=fit.f,replace=T),]))
	
	MN=genotype(parents[,2])+genotype(parents[,7])
	MT=parents[,8]
	GB.spread=0
	GB.u=(parents[,1]+parents[,6])/2
	GB=GB.u
	SEX=c(rep(1,nM),rep(0,nF))
	FIT=(1+ZMS*GB*SEX)*(1-abs(MN/2-MT)*MNS)
	pop=cbind(GB,MN,MT,SEX,FIT)

	return(pop)
}

parent=function(pop,i.gen){
	pop.m=pop[which(pop[,4]==1),]
	pop.f=pop[which(pop[,4]==0),]
	fit.m=pop.m[,5]/sum(pop.m[,5])
	fit.f=pop.f[,5]/sum(pop.f[,5])
	
	parents=rbind(cbind(pop.m[sample(1:nM,size=nM,prob=fit.m,replace=T),],pop.f[sample(1:nF,size=nM,prob=fit.f,replace=T),]),cbind(pop.m[sample(1:nM,size=nF,prob=fit.m,replace=T),],pop.f[sample(1:nF,size=nF,prob=fit.f,replace=T),]))
	
	return(parents)
}

evolve=function(pop0,FM0,PZ,MNS,ZMS){
	pop=pop0
	tSS=list()
	#tSS[['g0']]=c(mean(pop[,1]),mean(pop[,3]))
	pop=evolve_first(pop=pop,FM0=FM0,PZ=PZ,MNS=MNS,ZMS=ZMS) #First generation
	#tSS[['g1']]=c(mean(pop[,1]),mean(pop[,3]))
	i.gen=2
	while (i.gen <= nGEN){
		parents=parent(pop,i.gen)
		MN=genotype(parents[,2])+genotype(parents[,7])
		MT=parents[,8]
		GB.spread=sqrt(1/(60*(i.gen-1)))
		GB.u=(parents[,1]+parents[,6])/2
		GB=GB.u+rnorm(length(GB.u),mean=0,sd=GB.spread*GB.u*(1-GB.u))
		GB[which(GB<0)]=0 ; GB[which(1<GB)]=1
		SEX=c(rep(1,nM),rep(0,nF))
		FIT=(1+ZMS*GB*SEX)*(1-abs(MN/2-MT)*MNS)
		pop=cbind(GB,MN,MT,SEX,FIT)
		#if(i.gen<10 | i.gen%%10==0){
		#	tSS[[paste('g',i.gen,sep='')]]=c(mean(pop[,1]),mean(pop[,3]))	
		#}
		i.gen=i.gen + 1
	}
	#tSS[[paste('g',i.gen-1,sep='')]] = c(mean(pop[,1]),mean(pop[,2])/2,mean(pop[,3]))	
	#stats = unlist(tSS)
	stats = c(mean(pop[,1]),mean(pop[,2])/2,mean(pop[,3]))
	return(stats)
}

simulate=function(){	
	# Generating initial populations
	FM0=runif(1, min=0.01, max=max.FM0)
	PZ=runif(1, min=0.01, max=max.PZ)
	MNS=runif(1, min=0, max=max.MNS)
	ZMS=runif(1, min=0, max=max.ZMS)
	
	prop.zm=round(PS*PZ*FM0,0)
	prop.tm=round(PS*(1-PZ)*0.2,0)
	prop.zf=round(PS*PZ*(1-FM0),0)
	prop.tf=round(PS*(1-PZ)*0.8,0)
	
	# Individual coded as ('GB','MN','MT','SEX','FIT')
	zm=c(1,2,1,1,(1+ZMS*1*1)*(1))
	tm=c(0,0,0,1,(1+ZMS*1*0)*(1))
	zf=c(1,2,1,0,(1+ZMS*0*1)*(1))
	tf=c(0,0,0,0,(1+ZMS*0*0)*(1))
	
	pop=rbind(rbind(t(matrix(nrow=5,ncol=prop.zm,zm)), t(matrix(nrow=5,ncol=prop.tm,tm))),
	rbind(t(matrix(nrow=5,ncol=prop.zf,zf)),t(matrix(nrow=5,ncol=prop.tf,tf))))

	colnames(pop)=c('GB','MN','MT','SEX','FIT')
	stats = evolve(pop0=pop,FM0=FM0,PZ=PZ,MNS=MNS,ZMS=ZMS)
	tRESULT=c(FM0,PZ,MNS,ZMS,stats)
	return(tRESULT)
	}

main=function(){	
	tSUMMARY=t(replicate(n=nSET, expr=simulate(), simplify=T))
	print(paste('THREAD ',THR_NUM,': ',nSET,' population evolved',sep=''))
	#colnames(tSUMMARY)=c('FM0','PZ','MNS','ZMS',paste(c('GB','MT'),rep(c(seq(0,10),seq(20,151,10)),each=2),sep='.'),'MN')
	colnames(tSUMMARY)=c('FM0','PZ','MNS','ZMS','GB','MT','MN')
	save(tSUMMARY,file=paste('Rdata/set',THR_NUM,'.rda',sep=''))
}

####################################################################################
######################################MAIN RUN######################################
####################################################################################
args=commandArgs(trailingOnly=TRUE)
THR_NUM <<- as.double(args[1])
nSET <<- as.double(args[2])
nREP <<- as.double(args[3])
nM <<- as.double(args[4])
nF <<- as.double(args[5])
nGEN <<- as.integer(args[6])

max.FM0 <<- 0.99
max.PZ  <<- 0.99
max.MNS <<- 0.2
max.ZMS <<- 50
PS=nM+nF

main()



