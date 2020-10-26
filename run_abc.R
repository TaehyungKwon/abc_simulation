####################################################################################
############################MADE BY TAEHYUNG KWON Nov 11 2019#######################
#########################UPDATED BY TAEHYUNG KWON Mar 25 2020#######################
####################################################################################

## 1.Loading simulation set 
retrieve_data=function(){
	set=1
	tBOTH=list()
	tZS=list()
	tMN=list()
	while (set <= TOTAL.SET){ 
		infile1=paste('Rdata/a.set',set,'.rda',sep='')
		infile2=paste('Rdata/n.set',set,'.rda',sep='')
		infile3=paste('Rdata/n2.set',set,'.rda',sep='')
		if (file.exists(infile1)){load(infile1);tBOTH[[set]]=tSUMMARY}
		if (file.exists(infile2)){load(infile2);tZS[[set]]=tSUMMARY}
		if (file.exists(infile3)){load(infile3);tMN[[set]]=tSUMMARY}
		set=set + 1
	}
	print('all simulation set loaded')
	
	mZS<<-do.call('rbind',tZS)
	mMN<<-do.call('rbind',tMN)
	mBOTH<<-do.call('rbind',tBOTH)
	
	obs.ss<<-c(GB_obs,MT_obs,Y_obs)
	both.params<<-mBOTH[,1:5]
	both.sim.ss<<-mBOTH[,c(6,8,9)]
	zs.params<<-mZS[,1:5]
	zs.sim.ss<<-mZS[,c(6,8,9)]
	mn.params<<-mMN[,1:5]
	mn.sim.ss<<-mMN[,c(6,8,9)]
}

## 2.Running ABC
run_abc=function(){	
	
	if (file.exists(paste('Rdata/abc.',nSIM,'.rda',sep=''))==FALSE){
		tABC = list()
		tABC[['zs']]=list()
		tABC[['mn']]=list()
		tABC[['both']]=list()
		for (TOL in vTOLS){
			start.time=Sys.time()
			print(TOL)
			ctol=as.character(TOL)
			tABC[[1]][[ctol]] = abc(target=obs.ss, param=zs.params, sumstat=zs.sim.ss,tol=TOL,method="rejection")
			tABC[[2]][[ctol]] = abc(target=obs.ss, param=mn.params, sumstat=mn.sim.ss,tol=TOL,method="rejection")
			tABC[[3]][[ctol]] = abc(target=obs.ss, param=both.params, sumstat=both.sim.ss,tol=TOL,method="rejection")
			end.time=Sys.time()			
			print(paste('ABC performed: ',end.time-start.time,sep=''))
		}
		save(tABC, file=paste('Rdata/abc.',nSIM,'.rda',sep=''))
	} else {
		print(paste('ABC skipped: ABC result exists!'))
	}
	if (file.exists(paste('Rdata/abc.',nSIM,'.rda',sep=''))==T){
		load(file=paste('Rdata/abc.',nSIM,'.rda',sep=''))
	}
	
	tABC<<-tABC
}

## 3. Checking Goodness of Fit
run_gfit=function(){	 
	
	if (file.exists(paste('Rdata/gfit.',nSIM,'.rda',sep=''))==FALSE){
		tGFIT = list()
		tGFIT[['zs']]=list()
		tGFIT[['mn']]=list()
		tGFIT[['both']]=list()
		TOL=vTOLS[length(vTOLS)]
		start.time=Sys.time()
		ctol=as.character(TOL)
		tGFIT[[1]] = gfit(target=obs.ss, sumstat=zs.sim.ss, statistic=mean,tol=vTOLS[length(vTOLS)],nb.replicate=1000)
		tGFIT[[2]] = gfit(target=obs.ss, sumstat=mn.sim.ss, statistic=mean,tol=vTOLS[length(vTOLS)],nb.replicate=1000)
		tGFIT[[3]] = gfit(target=obs.ss, sumstat=both.sim.ss, statistic=mean,tol=vTOLS[length(vTOLS)],nb.replicate=1000)
		end.time=Sys.time()			
		print(paste('GFIT performed: ',end.time-start.time,sep=''))
		save(tGFIT, file=paste('Rdata/gfit.',nSIM,'.rda',sep=''))
	} else {
		print(paste('GFIT skipped: GFIT result exists!'))		
	}
	if (file.exists(paste('Rdata/gfit.',nSIM,'.rda',sep=''))==T){
		load(file=paste('Rdata/gfit.',nSIM,'.rda',sep=''))
	}
	
	tGFIT<<-tGFIT
	
}

## 4. Summarize the result
summarize=function(){
	library(ggplot2)
	library(cowplot)
	#library(egg)
	#library(reshape2)
	
	zs.abc=tABC[[1]][[length(vTOLS)]]
	mn.abc=tABC[[2]][[length(vTOLS)]]
	both.abc=tABC[[3]][[length(vTOLS)]]
		
	## Goodness of fit table
	vTYPE=c('male-biased zebu selection','mito-nuclear selection','both selections')
	for (i in 1:3){
		pdf(paste('results/gfit.',i,'.summary.pdf',sep=''))
		plot(tGFIT[[i]], main=paste(vTYPE[i],summary(tGFIT[[i]])$pvalue,sep=''))
		dev.off()
	}
		
	## Distance table
	dist.vec=c(tABC[[1]][[length(vTOLS)]]$dist,tABC[[2]][[length(vTOLS)]]$dist,tABC[[3]][[length(vTOLS)]]$dist)
	acc.vec=rep(0,nrow(mZS)+nrow(mMN)+nrow(mBOTH))
	for (i in 1:length(vTOLS)){
		acc.vec[which(tABC[[1]][[i]]$region)]=acc.vec[which(tABC[[1]][[i]]$region)]+1
		acc.vec[nrow(mZS)+which(tABC[[2]][[i]]$region)]=acc.vec[nrow(mZS)+which(tABC[[2]][[i]]$region)]+1
		acc.vec[nrow(mZS)+nrow(mMN)+which(tABC[[3]][[i]]$region)]=acc.vec[nrow(mZS)+nrow(mMN)+which(tABC[[3]][[i]]$region)]+1
	}
	type.vec=c(rep('zs',nrow(mZS)),rep('mn',nrow(mMN)),rep('both',nrow(mBOTH)))
	temp.df=data.frame(as.double(dist.vec),as.integer(acc.vec),as.factor(type.vec))
	colnames(temp.df)=c('dist','acc','type')
	
	df.dist=data.frame(c(rep(c('zs','mn','both'),each=4)),
	c(mean(temp.df[which(temp.df$type=='zs' & temp.df$acc>=0),1]),mean(temp.df[which(temp.df$type=='zs' & temp.df$acc>0),1]),mean(temp.df[which(temp.df$type=='zs' & temp.df$acc>1),1]),mean(temp.df[which(temp.df$type=='zs' & temp.df$acc>2),1]),
	mean(temp.df[which(temp.df$type=='mn' & temp.df$acc>=0),1]),mean(temp.df[which(temp.df$type=='mn' & temp.df$acc>0),1]),mean(temp.df[which(temp.df$type=='mn' & temp.df$acc>1),1]),mean(temp.df[which(temp.df$type=='mn' & temp.df$acc>2),1]),
	mean(temp.df[which(temp.df$type=='both' & temp.df$acc>=0),1]),mean(temp.df[which(temp.df$type=='both' & temp.df$acc>0),1]),mean(temp.df[which(temp.df$type=='both' & temp.df$acc>1),1]),mean(temp.df[which(temp.df$type=='both' & temp.df$acc>2),1])),
	c(sd(temp.df[which(temp.df$type=='zs' & temp.df$acc>=0),1]),sd(temp.df[which(temp.df$type=='zs' & temp.df$acc>0),1]),sd(temp.df[which(temp.df$type=='zs' & temp.df$acc>1),1]),sd(temp.df[which(temp.df$type=='zs' & temp.df$acc>2),1]),
	sd(temp.df[which(temp.df$type=='mn' & temp.df$acc>=0),1]),sd(temp.df[which(temp.df$type=='mn' & temp.df$acc>0),1]),sd(temp.df[which(temp.df$type=='mn' & temp.df$acc>1),1]),sd(temp.df[which(temp.df$type=='mn' & temp.df$acc>2),1]),
	sd(temp.df[which(temp.df$type=='both' & temp.df$acc>=0),1]),sd(temp.df[which(temp.df$type=='both' & temp.df$acc>0),1]),sd(temp.df[which(temp.df$type=='both' & temp.df$acc>1),1]),sd(temp.df[which(temp.df$type=='both' & temp.df$acc>2),1])))
	colnames(df.dist)=c('type','dist.mean','dist.stdev')
	
	write.table(file=paste('results/dist_summary.txt',sep=''),df.dist,quote=F,col.names=T,row.names=F,sep='\t')
		
	## Sumstat plot
	tHYP=list('zs'=mZS,'mn'=mMN,'both'=mBOTH)
	type.vec=c('zs','mn','both')
	tSS=list()
	for (i in c(1,3)){
		mat=cbind(tHYP[[i]][,6:9],rep(type.vec[i],nrow(tHYP[[i]])))
		acc.vec=rep('-',nrow(mat))
		for (j in 1:3){
			acc.vec[which(tABC[[i]][[j]]$region)]=names(tABC[[i]])[j]
		}
		mat.acc=cbind(mat[which(acc.vec=='0.001'),],acc.vec[which(acc.vec=='0.001')])
		tSS[[i]]=mat.acc
	}
	
	df.ss=data.frame(do.call(rbind,tSS))
	
	fig.ss=ggplot(df.ss)+
		geom_point(x=GB_obs,y=MT_obs,color='firebrick1', size=1, shape=18)+
		stat_density_2d(geom='polygon',aes(x=as.double(as.character(GB)),y=as.double(as.character(MT)),group=V5,fill=V5,alpha=..level../sum(..level..)),color='black',size=0.1)+
		scale_fill_manual(values=c('zs'='black','mn'='firebrick1','both'='deepskyblue1'))+
		labs(x='GB',y='MT')+
		theme(panel.background=element_rect(fill='white', colour='black',size=0.5),
			text=element_text(size=6,hjust=.5,vjust=.5),
			line=element_line(size=0.5),
			panel.border=element_rect(color='black',linetype='solid',fill=NA),
			panel.grid.major=element_line(color='grey',size=0.5),
			panel.grid.minor=element_blank(),
			axis.text.y=element_text(size=6,hjust=.5,vjust=.5),
			axis.text.x=element_text(size=6,hjust=.5,vjust=.5,angle=90),	
			plot.title=element_text(size=6,face='bold.italic'),
			plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
			axis.line=element_blank(),
			axis.ticks=element_line(size=0.25),
			strip.text=element_text(colour="white",size=8,margin=margin(0,0,0,0,'mm'),face="bold"),
			strip.background=element_rect(fill="#1E403C",color="black",linetype='solid',size=0.5),
			legend.background=element_rect(fill='white',size=0.5),
			legend.margin=margin(unit(0.1,'cm')),
			legend.key=element_rect(color='black'),
			legend.key.size=unit(0.2,'lines'),
			legend.position='right')
	
	ggsave(paste('results/sumstat.pdf',sep=''),fig.ss,units='mm',height=50,width=60,dpi=600)
	
	## Parameter inference: summary
	tSUM=list()
	tSUM[['zs']]=list()
	tSUM[['mn']]=list()
	tSUM[['both']]=list()
	for (tol in vTOLS){
		ctol=as.character(tol)
		tSUM[[1]][[ctol]]=summary(tABC[[1]][[ctol]],intvl=0.9)[c(2,3,6),]
		tSUM[[2]][[ctol]]=summary(tABC[[2]][[ctol]],intvl=0.9)[c(2,3,6),]
		tSUM[[3]][[ctol]]=summary(tABC[[3]][[ctol]],intvl=0.9)[c(2,3,6),]
	}
	tSUM<<-tSUM
	models=rep(c('zs','mn','both'),each=3*length(vTOLS))
	tolerance=rep(vTOLS,each=3,3)
	mSUM=cbind(models,tolerance,rbind(do.call(rbind,tSUM[[1]]),do.call(rbind,tSUM[[2]]),do.call(rbind,tSUM[[3]])))
	
	write.table(file=paste('results/params_summary.txt',sep=''),mSUM,quote=F,col.names=T,row.names=T,sep='\t')
	
	## Parameter inference: posterior and prior
	max_params=c(max_Fzm,max_Fzf,max_MF,max_Smn,max_Szs)
	df.both=data.frame(mBOTH)[both.abc$region,]
	
	temp.vec=rep('p',nrow(mBOTH))
	temp.vec[which(both.abc$region)]='a'
	df.acc=data.frame(mBOTH,temp.vec)
	
	fig.blank=ggplot()+geom_blank(aes(1,1))+
		theme(plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.margin=unit(c(0.01,0.01,0.01,0.01),"cm"),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank(),
		legend.position='none')
	
	tFIG=vector('list',length=25)
	for (i in 1:5){
		for (j in 1:5){
			index = 5*j-i+1
			if (i<j){
				temp.df=df.both[,c(i,j)]
				colnames(temp.df)=c('V1','V2')
				fig.jp=ggplot(temp.df,aes(x=V1, y=V2))+
					stat_density_2d(geom='polygon',aes(alpha=after_stat(level)))+
					geom_vline(xintercept=c(tSUM[[3]][[length(vTOLS)]][2,i]),color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					geom_hline(yintercept=c(tSUM[[3]][[length(vTOLS)]][2,j]),color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					scale_x_continuous(expand=c(0.005,0),limits=c(0,max_params[i]))+
					scale_y_continuous(expand=c(0.005,0),limits=c(0,max_params[j]))+
					labs(x='',y='',title='')+
					theme(text=element_text(size=6,hjust=.5,vjust=.5),
						line=element_line(size=0.5),
						panel.background=element_rect(fill='white', colour='black'),
						panel.border=element_rect(color='black',linetype='solid',fill=NA),
						panel.grid.minor=element_blank(),
						panel.grid.major=element_line(size=0.1,color='grey'),
						axis.text=element_blank(),
						axis.title=element_blank(),
						axis.ticks=element_line(size=0.25),
						plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"),
						axis.line=element_blank(),
						legend.background=element_rect(fill='white',size=0.5),
						legend.margin=margin(unit(0.1,'cm')),
						legend.key=element_rect(color='black'),
						legend.key.size=unit(1,'lines'),
						legend.position='none')
						
				tFIG[[index]]=fig.jp
				print(index)
				
			} else if (i==j) {	
				temp.acc=df.acc[,c(i,10)]
				colnames(temp.acc)=c('V1','V2')				
				fig.pp=ggplot(temp.acc,aes(x=V1,group=V2,fill=V2))+
					geom_vline(xintercept=tSUM[[3]][[length(vTOLS)]][2,i],color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					geom_vline(xintercept=c(tSUM[[3]][[length(vTOLS)]][1,i],tSUM[[3]][[length(vTOLS)]][3,i]),color='firebrick1',linetype='longdash',alpha=1,size=0.25)+
					geom_density(alpha=0.5,color=NA)+
					scale_x_continuous(expand=c(0.005,0),limits=c(0,max_params[i]))+
					scale_y_continuous(expand=c(0.005,0))+
					#scale_y_continuous(expand=c(0.005,0),limits=c(0,max_params[j]))+
					scale_fill_manual(values=c('p'='grey','a'='deepskyblue'))+
					labs(x='',y='',title='')+
					theme(text=element_text(size=6,hjust=.5,vjust=.5),
						line=element_line(size=0.5),
						panel.background=element_rect(fill='white', colour='black'),
						panel.border=element_rect(color='black',linetype='solid',fill=NA),
						panel.grid.minor=element_blank(),
						panel.grid.major=element_blank(),
						axis.text=element_blank(),
						axis.ticks=element_line(size=0.25),
						axis.title=element_blank(),
						plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"),
						axis.line=element_blank(),
						legend.background=element_rect(fill='white',size=0.5),
						legend.margin=margin(unit(0.1,'cm')),
						legend.key=element_rect(color='black'),
						legend.key.size=unit(1,'lines'),
						legend.position='none')
			
				tFIG[[index]]=fig.pp
			} else if (i>j){
				tFIG[[index]]=fig.blank
			}
		}
	}

	posteriorplot=plot_grid(plotlist=tFIG,ncol=5,nrow=5)
	save_plot(paste('results/ppdist.pdf',sep=''),posteriorplot,units='mm',base_height=150,base_width=150,dpi=1200)
		
	## Supplementary: posterior and prior for null distribution
	max_params=c(max_Fzm,max_Fzf,max_MF,max_Smn,max_Szs)
	df.zs=data.frame(mZS)[zs.abc$region,]
	
	temp.vec=rep('p',nrow(mZS))
	temp.vec[which(zs.abc$region)]='a'
	df.acc=data.frame(mZS,temp.vec)
	
	fig.blank=ggplot()+geom_blank(aes(1,1))+
		theme(plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.margin=unit(c(0.01,0.01,0.01,0.01),"cm"),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank(),
		legend.position='none')
	
	tFIG=vector('list',length=25)
	for (i in 1:5){
		for (j in 1:5){
			index = 5*j-i+1
			if (i<j){
				temp.df=df.zs[,c(i,j)]
				colnames(temp.df)=c('V1','V2')
				fig.jp=ggplot(temp.df,aes(x=V1, y=V2))+
					stat_density_2d(geom='polygon',aes(alpha=after_stat(level)))+
					geom_vline(xintercept=c(tSUM[[1]][[length(vTOLS)]][2,i]),color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					geom_hline(yintercept=c(tSUM[[1]][[length(vTOLS)]][2,j]),color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					scale_x_continuous(expand=c(0.005,0),limits=c(0,max_params[i]))+
					scale_y_continuous(expand=c(0.005,0),limits=c(0,max_params[j]))+
					labs(x='',y='',title='')+
					theme(text=element_text(size=6,hjust=.5,vjust=.5),
						line=element_line(size=0.5),
						panel.background=element_rect(fill='white', colour='black'),
						panel.border=element_rect(color='black',linetype='solid',fill=NA),
						panel.grid.minor=element_blank(),
						panel.grid.major=element_line(size=0.1,color='grey'),
						axis.text=element_blank(),
						axis.title=element_blank(),
						axis.ticks=element_line(size=0.25),
						plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"),
						axis.line=element_blank(),
						legend.background=element_rect(fill='white',size=0.5),
						legend.margin=margin(unit(0.1,'cm')),
						legend.key=element_rect(color='black'),
						legend.key.size=unit(1,'lines'),
						legend.position='none')
						
				tFIG[[index]]=fig.jp
				print(index)
				
			} else if (i==j) {	
				temp.acc=df.acc[,c(i,10)]
				colnames(temp.acc)=c('V1','V2')				
				fig.pp=ggplot(temp.acc,aes(x=V1,group=V2,fill=V2))+
					geom_vline(xintercept=tSUM[[1]][[length(vTOLS)]][2,i],color='firebrick1',linetype='solid',alpha=1,size=0.25)+
					geom_vline(xintercept=c(tSUM[[1]][[length(vTOLS)]][1,i],tSUM[[1]][[length(vTOLS)]][3,i]),color='firebrick1',linetype='longdash',alpha=1,size=0.25)+
					geom_density(alpha=0.5,color=NA)+
					scale_x_continuous(expand=c(0.005,0),limits=c(0,max_params[i]))+
					scale_y_continuous(expand=c(0.005,0))+
					#scale_y_continuous(expand=c(0.005,0),limits=c(0,max_params[j]))+
					scale_fill_manual(values=c('p'='grey','a'='deepskyblue'))+
					labs(x='',y='',title='')+
					theme(text=element_text(size=6,hjust=.5,vjust=.5),
						line=element_line(size=0.5),
						panel.background=element_rect(fill='white', colour='black'),
						panel.border=element_rect(color='black',linetype='solid',fill=NA),
						panel.grid.minor=element_blank(),
						panel.grid.major=element_blank(),
						axis.text=element_blank(),
						axis.ticks=element_line(size=0.25),
						axis.title=element_blank(),
						plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"),
						axis.line=element_blank(),
						legend.background=element_rect(fill='white',size=0.5),
						legend.margin=margin(unit(0.1,'cm')),
						legend.key=element_rect(color='black'),
						legend.key.size=unit(1,'lines'),
						legend.position='none')
			
				tFIG[[index]]=fig.pp
			} else if (i>j){
				tFIG[[index]]=fig.blank
			}
		}
	}

	posteriorplot=plot_grid(plotlist=tFIG,ncol=5,nrow=5)

	save_plot(paste('results/ppdist_zs.pdf',sep=''),posteriorplot,units='mm',base_height=150,base_width=150,dpi=1200)
	
	## Supplementary: GB and MN
	fig.ssmn=ggplot(df.ss)+
		stat_density_2d(geom='polygon',aes(x=as.double(as.character(GB)),y=as.double(as.character(MN)),group=V5,fill=V5,alpha=..level../sum(..level..)),color='black',size=0.1)+
		scale_fill_manual(values=c('zs'='black','mn'='firebrick1','both'='deepskyblue1'))+
		labs(x='GB',y='MN')+
		theme(panel.background=element_rect(fill='white', colour='black',size=0.5),
			text=element_text(size=6,hjust=.5,vjust=.5),
			line=element_line(size=0.5),
			panel.border=element_rect(color='black',linetype='solid',fill=NA),
			panel.grid.major=element_line(color='grey',size=0.5),
			panel.grid.minor=element_blank(),
			axis.text.y=element_text(size=6,hjust=.5,vjust=.5),
			axis.text.x=element_text(size=6,hjust=.5,vjust=.5,angle=90),	
			plot.title=element_text(size=6,face='bold.italic'),
			plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
			axis.line=element_blank(),
			axis.ticks=element_line(size=0.25),
			strip.text=element_text(colour="white",size=8,margin=margin(0,0,0,0,'mm'),face="bold"),
			strip.background=element_rect(fill="#1E403C",color="black",linetype='solid',size=0.5),
			legend.background=element_rect(fill='white',size=0.5),
			legend.margin=margin(unit(0.1,'cm')),
			legend.key=element_rect(color='black'),
			legend.key.size=unit(0.2,'lines'),
			legend.position='right')
	
	ggsave(paste('results/sumstat.mn_gb.pdf',sep=''),fig.ssmn,units='mm',height=50,width=60,dpi=600)
	
}

temp_call=function(){

	##GLOBAL PARAMETERS
	nSIM=10000000
	nTHREAD=100
	nSET=1000

	##SIMULATION PARAMETERS
	nPOP = 10000	#Population size at nGEN
	nGEN = 110		#Number of generation, approximated by ALDER
	max_Fzm = 1
	max_Fzf = 0.5
	max_MF = 0.5
	max_Smn = 0.2
	max_Szs = 100

	##ABC PARAMETERS
	GB_obs = 0.755
	MT_obs = 0
	Y_obs = 1
	TOL = 0.001		#Tolerence of ABC


	##ANNOUNCE 
	#GLOBAL PARAMETERS
	nSIM<<-nSIM
	nTHREAD<<-nTHREAD
	nSET<<-nSET

	#SIMULATION PARAMETERS
	nPOP<<-nPOP
	nGEN<<-nGEN
	max_Fzm<<-max_Fzm
	max_Fzf<<-max_Fzf
	max_MF<<-max_MF
	max_Smn<<-max_Smn
	max_Szs<<-max_Szs
		
	#ABC PARAMETERS
	GB_obs<<-GB_obs
	MT_obs<<-MT_obs
	Y_obs<<-Y_obs
	TOL<<-TOL
	TOTAL.SET<<-nSIM/nSET
}
#################################################################################
## ARGUMENTS
args=commandArgs(trailingOnly=TRUE)
TOTAL.SET<<-as.integer(args[1])
nSIM<<-as.integer(args[2])
GB_obs<<-as.double(args[3])
MT_obs<<-as.double(args[4])
Y_obs<<-as.double(args[5])
TOL<<-as.double(args[6])
max_Fzm<<-as.double(args[7])
max_Fzf<<-as.double(args[8])
max_MF<<-as.double(args[9])
max_Smn<<-as.double(args[10])
max_Szs<<-as.double(args[11])
vTOLS<<-c(0.01,0.005,0.001)

## FUNCTIONS
## Requisite library: abc, ggplot2, cowplot
library(abc)
retrieve_data()
run_abc()
run_gfit()
summarize()
