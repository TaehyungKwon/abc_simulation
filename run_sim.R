rm(list=ls())
####################################################################################
###########################MADE BY TAEHYUNG KWON Nov 11 2019########################
########################CORRECTED BY TAEHYUNG KWON MFr 18 2020######################
####################################################################################
meiosis_haplo=function(geno){
	size=nrow(geno)
	hap.vec=c(geno[,1],geno[,2])
	select.hap=sample(c(0,size),size,0.5)
	geno.hap=hap.vec[1:size+select.hap]
	return(geno.hap)
}

meiosis_recomb=function(geno,nsplit){
	size=nrow(geno)
	hap.u=(geno[,1]+geno[,2])/2
	geno.hap=rnorm(n=size,mean=hap.u,sd=sqrt(hap.u*(1-hap.u)/nsplit))
	geno.hap[which(geno.hap<0)]=0 ; geno.hap[which(1<geno.hap)]=1
	return(geno.hap)
}

evolve.gen1=function(pop,MF,Smn,Szs){
	## When all parents have homozygote GB & MN
	pop.m=pop[which(pop[,7]==1),]
	pop.f=pop[which(pop[,7]==0),]
	fit.m=pop.m[,8]
	fit.f=pop.f[,8]
	parents=cbind(pop.m[sample(1:nrow(pop.m),size=nPOP*dN[2],prob=fit.m,replace=T),],pop.f[sample(1:nrow(pop.f),size=nPOP*dN[2],prob=fit.f,replace=T),])
	
	nM=sum(sample(c(0,1),size=nPOP*dN[2],replace=T,prob=c(1-MF,MF)))
	nF=nPOP*dN[2]-nM
	GB.1=parents[,1]
	GB.2=parents[,9]
	GB.u=(GB.1+GB.2)/2
	MN.1=parents[,3]
	MN.2=parents[,11]
	MN.u=(MN.1+MN.2)/2
	MT=parents[,13]
	Y=parents[,6]
	SEX=c(rep(1,nM),rep(0,nF))
	FIT=(1+Szs*(GB.u*Y)*SEX)*(1-abs(MN.u-MT)*Smn)
	pop=cbind(GB.1,GB.2,MN.1,MN.2,MT,Y,SEX,FIT)
	return(pop)
}

evolve=function(pop0,MF,Smn,Szs){
	pop=pop0
	pop=evolve.gen1(pop=pop,MF=MF,Smn=Smn,Szs=Szs) 
	i.gen=2
	
	while (i.gen <= nGEN){
		## Parent sampling
		pop.m=pop[which(pop[,7]==1),]
		pop.f=pop[which(pop[,7]==0),]
		fit.m=pop.m[,8]
		fit.f=pop.f[,8]
		parents=cbind(pop.m[sample(1:nrow(pop.m),size=nPOP*dN[i.gen+1],prob=fit.m,replace=T),],pop.f[sample(1:nrow(pop.f),size=nPOP*dN[i.gen+1],prob=fit.f,replace=T),])
	
		## Meiosis and mating
		nM=sum(sample(c(0,1),size=nPOP*dN[i.gen+1],replace=T,prob=c(1-MF,MF)))
		nF=nPOP*dN[i.gen+1]-nM
		GB.1=meiosis_recomb(parents[,c(1,2)],25*(i.gen))
		GB.2=meiosis_recomb(parents[,c(9,10)],25*(i.gen))
		GB.u=(GB.1+GB.2)/2
		MN.1=meiosis_haplo(parents[,c(3,4)])
		MN.2=meiosis_haplo(parents[,c(11,12)])
		MN.u=(MN.1+MN.2)/2
		MT=parents[,13]
		Y=parents[,6]
		SEX=c(rep(1,nM),rep(0,nF))
		FIT=(1+Szs*(GB.u*Y)*SEX)*(1-abs(MN.u-MT)*Smn)
		pop=cbind(GB.1,GB.2,MN.1,MN.2,MT,Y,SEX,FIT)
		i.gen=i.gen+1
	}
	stats=c(mean(pop[,c(1,2)]),mean(pop[,c(3,4)]),mean(pop[,5]),mean(pop[,6]))
	return(stats)
}

simulate=function(){	
	## Sampling priors
	Fzm=round(runif(1, min=0, max=max_Fzm),4)
	Fzf=round(runif(1, min=0, max=max_Fzf),4)
	MF=round(runif(1, min=0.05, max=max_MF),4)
	Smn=round(runif(1, min=0, max=max_Smn),4)
	Szs=round(runif(1, min=0, max=max_Szs),2)
	
	## Generating each groups
	prop.zm=round(nPOP*dN[1]*MF*Fzm,0)
	prop.tm=round(nPOP*dN[1]*MF*(1-Fzm),0)
	prop.zf=round(nPOP*dN[1]*(1-MF)*Fzf,0)
	prop.tf=round(nPOP*dN[1]*(1-MF)*(1-Fzf),0)
	
	## Adjusting zero groups
	vPROP=c(prop.zm,prop.zf,prop.tm,prop.tf)	
	for (k in 1:4){
		if (vPROP[k]==0){
			vPROP[k]=1
		}
	}

	## Individuals are coded as ('GB.0','GB.1','MN.0','MN.1','MT','Y','SEX','FIT')
	zm=c(c(1,1),c(1,1),1,1,1,(1+Szs*(1*1)*1))
	tm=c(c(0,0),c(0,0),0,0,1,(1+Szs*(0*0)*1))
	zf=c(c(1,1),c(1,1),1,0,0,(1+Szs*(1*0)*0))
	tf=c(c(0,0),c(0,0),0,0,0,(1+Szs*(0*0)*0))
	
	pop=rbind(t(matrix(nrow=8,ncol=vPROP[1],zm)),t(matrix(nrow=8,ncol=vPROP[3],tm)),t(matrix(nrow=8,ncol=vPROP[2],zf)),t(matrix(nrow=8,ncol=vPROP[4],tf)))
	colnames(pop)=c('GB.0','GB.1','MN.0','MN.1','MT','Y','SEX','FIT')
	stats=evolve(pop0=pop,MF=MF,Smn=Smn,Szs=Szs)
	tRESULT=c(Fzm,Fzf,MF,Smn,Szs,stats)	
	return(tRESULT)
	}

popsize=function(){
	dN.N=1
	dN.67=1.031683
	dN.32=0.897195
	dN.16=0.681529
	dN.9=0.477029
	dN.4=0.323883
	dN.2=0.219904

	dN<<-c(rep(dN.N,nGEN-67),rep(dN.67,35),rep(dN.32,16),rep(dN.16,7),rep(dN.9,5),rep(dN.4,2),rep(dN.2,3))
}

## Execute simulation for alternative hypothesis
both=function(){
	popsize()
	outfile=paste('Rdata/a.set',THR_NUM,'.rda',sep='')
	if (!file.exists(outfile)){
		max_Smn<<-0
		max_Szs<<-as.double(args[9])
		tSUMMARY=t(replicate(n=nSET*0.1, expr=simulate(), simplify=T))
		colnames(tSUMMARY)=c('Fzm','Fzf','MF','Smn','Szs','GB','MN','MT','Y')
		save(tSUMMARY,file=outfile)
	}
}
## Execute simulation for null hypothesis
zs=function(){
	outfile=paste('Rdata/n.set',THR_NUM,'.rda',sep='')
	if (!file.exists(outfile)){
		max_Smn<<-0
		max_Szs<<-as.double(args[9])
		tSUMMARY=t(replicate(n=nSET*0.1, expr=simulate(), simplify=T))
		colnames(tSUMMARY)=c('Fzm','Fzf','MF','Smn','Szs','GB','MN','MT','Y')
		save(tSUMMARY,file=outfile)
	}
}

mn=function(){
	outfile=paste('Rdata/n2.set',THR_NUM,'.rda',sep='')
	if (!file.exists(outfile)){
		print('n2')
		max_Smn<<-as.double(args[8])
		max_Szs<<-0
		tSUMMARY=t(replicate(n=nSET*0.1, expr=simulate(), simplify=T))
		colnames(tSUMMARY)=c('Fzm','Fzf','MF','Smn','Szs','GB','MN','MT','Y')
		save(tSUMMARY,file=outfile)
	}
}


####################################################################################
######################################MAIN RUN######################################
####################################################################################
args=commandArgs(trailingOnly=TRUE)
THR_NUM<<-as.double(args[1])
nSET<<-as.double(args[2])
nGEN<<-as.integer(args[3])
nPOP<<-as.double(args[4])
max_Fzm<<-as.double(args[5])
max_Fzf<<-as.double(args[6])
max_MF<<-as.double(args[7])
max_Smn<<-as.double(args[8])
max_Szs<<-as.double(args[9])


both()
zs()
mn()

print(paste('THREAD ',THR_NUM,': FINISHED...',sep=''))