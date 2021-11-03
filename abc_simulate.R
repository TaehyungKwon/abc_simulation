### MADE BY TAEHYUNG KWON
in_popsize=function(){
	# popsize change was estimated using SMC++
	dNe.N=1
	dNe.67=1.031683
	dNe.32=0.897195
	dNe.16=0.681529
	dNe.9=0.477029
	dNe.4=0.323883
	dNe.2=0.219904
	dNe<<-c(rep(dNe.N,ngen-67),rep(dNe.67,35),rep(dNe.32,16),rep(dNe.16,7),rep(dNe.9,5),rep(dNe.4,2),rep(dNe.2,3))
}
read_params=function(){
	rm(list=ls())
	args=commandArgs(trailingOnly=TRUE)
	PRC<<-as.double(args[1])
	nREP<<-as.double(args[2])
	data.params=read.table(paste('Rdata/params.sim.txt',sep=''),header=F)
	for (row in 1:nrow(data.params)){
		assign(as.character(data.params[row,1]),data.params[row,2],envir=.GlobalEnv)
	}
	in_popsize()
}
in_meiosis_haplo=function(geno){
	size=nrow(geno)
	hap.vec=c(geno[,1],geno[,2])
	index.hap=sample(c(0,size),size,0.5)
	hap=hap.vec[1:size+index.hap]
	return(hap)
}
in_meiosis_recomb=function(geno,numfrag){
	size=nrow(geno)
	geno.u=(geno[,1]+geno[,2])/2
	hap=rnorm(n=size,mean=geno.u,sd=sqrt(geno.u*(1-geno.u)/numfrag))
	hap[which(hap<0)]=0;hap[which(1<hap)]=1
	return(hap)
}
in_initialize=function(pop,MF,Smn,Szs,Ne){
	## Gen 1/all parents have homozygote GB & MN
	pop.m=pop[which(pop[,7]==1),,drop=F]
	pop.f=pop[which(pop[,7]==0),,drop=F]
	## Population perished
	if ((nrow(pop.m)==0)|(nrow(pop.f)==0)) return(rep(-1,4))
	popsize=Ne*dNe[2]
	fit.m=pop.m[,8]
	fit.f=pop.f[,8]
	parents=cbind(pop.m[sample(1:nrow(pop.m),size=popsize,prob=fit.m,replace=T),],pop.f[sample(1:nrow(pop.f),size=popsize,prob=fit.f,replace=T),])
	GB.1=parents[,1]
	GB.2=parents[,9]
	GB.u=(GB.1+GB.2)/2
	MN.1=parents[,3]
	MN.2=parents[,11]
	MN.u=(MN.1+MN.2)/2
	MT=parents[,13]
	Y=parents[,6]
	SEX=sample(c(rep(0,round(popsize*(1-MF),0)),rep(1,round(popsize*MF,0))))
	FIT=(1-Smn*abs(MN.u-MT))*(1+Szs*(GB.u*Y)*SEX)
	pop=cbind(GB.1,GB.2,MN.1,MN.2,MT,Y,SEX,FIT)
	return(pop)
}
in_evolve=function(pop0,MF,Smn,Szs,Ne){
	pop=in_initialize(pop=pop0,MF=MF,Smn=Smn,Szs=Szs,Ne=Ne) 
	i.gen=2
	while (i.gen <= ngen){
		## Parent sampling
		pop.m=pop[which(pop[,7]==1),,drop=F]
		pop.f=pop[which(pop[,7]==0),,drop=F]
		popsize=Ne*dNe[i.gen+1]
		fit.m=pop.m[,8]
		fit.f=pop.f[,8]
		parents=cbind(pop.m[sample(1:nrow(pop.m),size=popsize,prob=fit.m,replace=T),],pop.f[sample(1:nrow(pop.f),size=popsize,prob=fit.f,replace=T),])
		## Meiosis and mating
		GB.1=in_meiosis_recomb(parents[,c(1,2)],25*(i.gen-1))
		GB.2=in_meiosis_recomb(parents[,c(9,10)],25*(i.gen-1))
		GB.u=(GB.1+GB.2)/2
		MN.1=in_meiosis_haplo(parents[,c(3,4)])
		MN.2=in_meiosis_haplo(parents[,c(11,12)])
		MN.u=(MN.1+MN.2)/2
		MT=parents[,13]
		Y=parents[,6]
		SEX=sample(c(rep(0,round(popsize*(1-MF),0)),rep(1,round(popsize*MF,0))))
		FIT=(1-Smn*abs(MN.u-MT))*(1+Szs*(GB.u*Y)*SEX)
		pop=cbind(GB.1,GB.2,MN.1,MN.2,MT,Y,SEX,FIT)
		i.gen=i.gen+1
	}
	GB.vec=(pop[,1]+pop[,2])/2
	MN.vec=(pop[,3]+pop[,4])/2
	MT.vec=pop[,5]
	Y.vec=pop[,6]
	stats=c(mean(MN.vec),
		mean(GB.vec),
		mean(MT.vec),
		mean(Y.vec))
	return(stats)
}
simulate=function(Ne,max.Smn,max.Szs,model){
	## Sampling priors
	Fzm=round(runif(1, min=0, max=max_Fzm),4)
	Fzf=round(runif(1, min=0, max=max_Fzf),4)
	MF=round(runif(1, min=min_MF, max=max_MF),4)
	Smn=round(runif(1, min=0, max=max.Smn),4)
	Szs=round(runif(1, min=0, max=max.Szs),2)
	## Generating each groups
	size.zm=round(Ne*dNe[1]*MF*Fzm,0)
	size.tm=round(Ne*dNe[1]*MF*(1-Fzm),0)
	size.zf=round(Ne*dNe[1]*(1-MF)*Fzf,0)
	size.tf=round(Ne*dNe[1]*(1-MF)*(1-Fzf),0)
	## Adjusting zero groups
	if ((size.zm+size.tm==0)|(size.zf+size.tf==0)) return(c(Fzm,Fzf,MF,Smn,Szs,Ne,rep(-1,4)))
	## Individuals have following properties: 'GB.0','GB.1','MN.0','MN.1','MT','Y','SEX','FIT'
	## FIT=(1–Smn*|MN – MT|)*(1+Szs*((GB*Y)*SEX))
	zm=c(c(1,1),c(1,1),1,1,1,(1+Szs*(1*1)*1))
	tm=c(c(0,0),c(0,0),0,0,1,(1+Szs*(0*0)*1))
	zf=c(c(1,1),c(1,1),1,0,0,(1+Szs*(1*0)*0))
	tf=c(c(0,0),c(0,0),0,0,0,(1+Szs*(0*0)*0))
	pop=rbind(t(matrix(nrow=8,ncol=size.zm,zm)),t(matrix(nrow=8,ncol=size.tm,tm)),t(matrix(nrow=8,ncol=size.zf,zf)),t(matrix(nrow=8,ncol=size.tf,tf)))
	colnames(pop)=c('GB.0','GB.1','MN.0','MN.1','MT','Y','SEX','FIT')
	stats=in_evolve(pop0=pop,MF=MF,Smn=Smn,Szs=Szs,Ne=Ne)
	tRESULT=c(Fzm,Fzf,MF,Smn,Szs,Ne,stats)
	return(tRESULT)
}
model_run=function(){
	read_params()
	for (model in c('neutral','mnsel','zmsel','bothsel')){
		if (model=='neutral') {
			max.Smn=0; max.Szs=0
		} else if (model=='zmsel') {
			max.Smn=0; max.Szs=max_Szs
		} else if (model=='mnsel') {
			max.Smn=max_Smn; max.Szs=0
		} else if (model=='bothsel') {
			max.Smn=max_Smn; max.Szs=max_Szs
		}
		for (Ne in unique(c(Ne_1,Ne_2,Ne_3))){
			outfile=paste('Rdata/',model,'.Ne',Ne,'.rep',PRC,'.rda',sep='')
			if (!file.exists(outfile)){
				tSUMMARY=t(replicate(n=nREP, expr=simulate(Ne=Ne,max.Smn=max.Smn,max.Szs=max.Szs,model=model), simplify=T))
				colnames(tSUMMARY)=c('Fzm','Fzf','MF','Smn','Szs','Ne','MN.mean','GB.mean','MT.mean','Y.mean')
				save(tSUMMARY,file=outfile)
			}
		}
	}
}
model_run()