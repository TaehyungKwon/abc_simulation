### MADE BY TAEHYUNG KWON
## Requisites: abc, ggplot2, cowplot
library(ggplot2)
library(cowplot)
library(abc)

deprecated_test=function(){
	require(abc.data)
	data(human)
	stat.voight
	TOTAL.REP<<-500
	nREP<<-10000
	nSIM<<-TOTAL.REP*nREP
	
	data.params.sim=read.table(paste('Rdata/params.sim.txt',sep=''),header=F)
	data.params.abc=read.table(paste('Rdata/params.abc.txt',sep=''),header=F)
	data.params=rbind(data.params.sim,data.params.abc)
	for (row in 1:nrow(data.params)){
		assign(as.character(data.params[row,1]),data.params[row,2],envir=.GlobalEnv)
	}
	Ne.vec<<-unique(as.character(c(Ne_1,Ne_2,Ne_3)))
	obs.sumstat<<-c(GB_mean,MT_mean,Y_mean)
	tol.vec<<-as.character(c(0.01,0.001))
	color.model<<-c('neutral'='grey60','mnsel'='dodgerblue2','zmsel'='firebrick3','bothsel'='purple2')

	model.vec<<-c('neutral','mnsel','zmsel','bothsel')
	tMATRIX=list()
	tPARAMS=list()
	tSUMSTAT=list()
	tSURVIVAL=c()
	for (Ne in Ne.vec){
		tMATRIX[[Ne]]=list()
		tPARAMS[[Ne]]=list()
		tSUMSTAT[[Ne]]=list()
		for (model in model.vec){
			temp.list=list()
			rep=1
			while (rep <= TOTAL.REP){ 
				infile=paste('Rdata/',model,'.Ne',Ne,'.rep',rep,'.rda',sep='')
				if (file.exists(infile)){load(infile);temp.list[[rep]]=tSUMMARY}
				rep=rep+1
			}
			temp.matrix=do.call('rbind',temp.list)
			temp.survived=temp.matrix[which(temp.matrix[,7]!=-1),]
			tMATRIX[[Ne]][[model]]=temp.survived
			tPARAMS[[Ne]][[model]]=temp.survived[,1:5]
			tSUMSTAT[[Ne]][[model]]=temp.survived[,8:ncol(temp.survived)]
		}
	}
	
	print(paste('INFO: ',TOTAL.REP,' replicates loaded',sep=''))
	tMATRIX<<-tMATRIX
	tPARAMS<<-tPARAMS
	tSUMSTAT<<-tSUMSTAT
	run_abc()
}
in_read_params=function(){
	rm(list=ls())
	args=commandArgs(trailingOnly=TRUE)
	TOTAL.REP<<-as.integer(args[1])
	nREP<<-as.integer(args[2])
	nSIM<<-TOTAL.REP*nREP
	data.params.sim=read.table(paste('Rdata/params.sim.txt',sep=''),header=F)
	data.params.abc=read.table(paste('Rdata/params.abc.txt',sep=''),header=F)
	data.params=rbind(data.params.sim,data.params.abc)
	for (row in 1:nrow(data.params)){
		assign(as.character(data.params[row,1]),data.params[row,2],envir=.GlobalEnv)
	}
	Ne.vec<<-unique(as.character(c(Ne_1,Ne_2,Ne_3)))
	obs.sumstat<<-c(GB_mean,MT_mean,Y_mean)
	tol.vec<<-as.character(c(0.01,0.001))
	color.model<<-c('neutral'='grey60','mnsel'='dodgerblue2','zmsel'='firebrick3','bothsel'='purple2')
}
in_visualize_population_survival=function(){
	i=1
	temp.list=list()
	for (Ne in Ne.vec){
		for (model in model.vec){
			temp.matrix=tMATRIX[[Ne]][[model]]
			temp.survived=temp.matrix[which(temp.matrix[,7]!=-1),]
			temp.perished=temp.matrix[which(temp.matrix[,7]==-1),]
			temp.list[[i]]=cbind(Ne=rep(Ne,2),model=rep(model,2),survival=c('survived','perished'),percentage=c(round(nrow(temp.survived)/nSIM*100,2),round((nSIM-nrow(temp.survived))/nSIM*100,2)))
			if (nrow(temp.survived)<100){
				print(paste('WARNING: population almost perished (<100): Ne=',Ne,'; model=',model,'; survived=',nrow(temp.survived),sep=''))
				next
			}
			i=i+1
		}
	}
	print('INFO: visualize population survival')
	df=data.frame(do.call(rbind,temp.list),stringsAsFactors=F)
	df$percentage=as.numeric(df$percentage)
	df$model=factor(df$model,levels=model.vec)
	df$Ne=factor(df$Ne,levels=Ne.vec)
	df$survival=factor(df$survival,levels=c('survived','perished'))
	df=df[which(df$survival=='perished'),]
	theme.common=theme(
		panel.background=element_rect(fill='white',color='black',size=0.5),
		text=element_text(size=7,hjust=.5,vjust=.5),
		line=element_line(size=0.5),
		panel.border=element_rect(color='black',linetype='solid',fill=NA),
		panel.grid.major=element_line(color='grey90',size=0.25),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=7,hjust=.5,vjust=.5),	
		plot.title=element_text(size=7,face='bold.italic'),
		plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		axis.line=element_blank(),
		axis.ticks=element_line(size=0.25),
		strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
		strip.background=element_rect(fill="grey90",color=NA,linetype='solid',size=0.5),
		legend.background=element_rect(fill='white',size=0.5),
		legend.margin=margin(unit(0.1,'cm')),
		legend.key=element_rect(color='black'),
		legend.key.size=unit(0.2,'lines'),
		legend.position='none')
	
	fig=ggplot(df,aes(x=Ne,y=percentage,fill=survival))+
		geom_bar(color=NA,position='stack',stat='identity',width=0.8)+
		scale_fill_manual(values=c('survived'='dodgerblue2','perished'='grey30'))+
		theme.common+
		facet_grid(.~model)
	ggsave(paste('out_figures/population_perished.pdf',sep=''),fig,units='mm',height=50,width=150,dpi=600)
}
in_read_data=function(){
	model.vec<<-c('neutral','mnsel','zmsel','bothsel')
	tMATRIX=list()
	tPARAMS=list()
	tSUMSTAT=list()
	tSURVIVAL=c()
	for (Ne in Ne.vec){
		tMATRIX[[Ne]]=list()
		tPARAMS[[Ne]]=list()
		tSUMSTAT[[Ne]]=list()
		for (model in model.vec){
			temp.list=list()
			rep=1
			while (rep <= TOTAL.REP){ 
				infile=paste('Rdata/',model,'.Ne',Ne,'.rep',rep,'.rda',sep='')
				if (file.exists(infile)){load(infile);temp.list[[rep]]=tSUMMARY}
				rep=rep+1
			}
			temp.matrix=do.call('rbind',temp.list)
			temp.survived=temp.matrix[which(temp.matrix[,7]!=-1),]
			tMATRIX[[Ne]][[model]]=temp.survived
			tPARAMS[[Ne]][[model]]=temp.survived[,1:5]
			tSUMSTAT[[Ne]][[model]]=temp.survived[,8:ncol(temp.survived)]
		}
	}
	print(paste('INFO: ',TOTAL.REP,' replicates loaded',sep=''))
	tMATRIX<<-tMATRIX
	tPARAMS<<-tPARAMS
	tSUMSTAT<<-tSUMSTAT
}
## Load parameters and data
load_data=function(){
	in_read_params()
	in_read_data()
	in_visualize_population_survival()
}
## Run ABC
run_abc=function(){	
	dir.create(file.path('Rdata/abc_result'), showWarnings=FALSE)
	if (!file.exists(paste('Rdata/abc_result/abc.',nSIM,'.rda',sep=''))){
		tABC=list()
		for (Ne in Ne.vec){
			tABC[[Ne]]=list()
			for (model in model.vec){
				tABC[[Ne]][[model]]=list()
				for (tol in tol.vec){
					tABC[[Ne]][[model]][[tol]]=abc(target=obs.sumstat, param=tPARAMS[[Ne]][[model]], sumstat=tSUMSTAT[[Ne]][[model]],tol=as.numeric(tol),method="rejection")
				}
			}
		}
		print(paste('INFO: perform ABC',sep=''))
		save(tABC, file=paste('Rdata/abc_result/abc.',nSIM,'.rda',sep=''))
	} else {
		print(paste('INFO: load previous ABC result'))
		load(file=paste('Rdata/abc_result/abc.',nSIM,'.rda',sep=''))
	}
	tABC<<-tABC
}
## Run model selection
run_model_selection=function(){
	if (!file.exists(paste('Rdata/abc_result/modsel.',nSIM,'.rda',sep=''))){
		tMODSEL=list()
		for (Ne in Ne.vec){
			print('INFO: perform MODSEL')
			tMODSEL[[Ne]]=list()
			model.concat=c()
			for (model in model.vec){
				model.concat=c(model.concat,rep(model,nrow(tSUMSTAT[[Ne]][[model]])))
			}
			sumstat.concat=do.call(rbind,tSUMSTAT[[Ne]])
			for (tol in tol.vec){
				tMODSEL[[Ne]][[tol]]=list()
				tMODSEL[[Ne]][[tol]][['cvmodsel']]=cv4postpr(model.concat, sumstat.concat, nval=20, tol=as.numeric(tol), method="rejection")
				tMODSEL[[Ne]][[tol]][['modsel']]=postpr(obs.sumstat, model.concat, sumstat.concat, tol=as.numeric(tol), method="rejection")
			}
		}
		save(tMODSEL, file=paste('Rdata/abc_result/modsel.',nSIM,'.rda',sep=''))
	} else {
		print(paste('INFO: load previous MODSEL result'))
		load(file=paste('Rdata/abc_result/modsel.',nSIM,'.rda',sep=''))
	}
	tMODSEL<<-tMODSEL
}
## Run goodness-of-fit test
run_goodness_of_fit=function(){
	bool.gfit=TRUE
	if (bool.gfit){
		if (!file.exists(paste('Rdata/abc_result/gfit.',nSIM,'.rda',sep=''))){
			tGFIT=list()
			start.time=Sys.time()
			for (Ne in Ne.vec){
				tGFIT[[Ne]]=list()
				for (model in model.vec){
					tGFIT[[Ne]][[model]]=list()
					tGFIT[[Ne]][[model]]=gfit(target=obs.sumstat,
						sumstat=tSUMSTAT[[Ne]][[model]], statistic=mean, nb.replicate=100)
				}
			}
			end.time=Sys.time()
			print(paste('INFO: perform Goodness-of-fit test',sep=''))
			save(tGFIT, file=paste('Rdata/abc_result/gfit.',nSIM,'.rda',sep=''))
		} else {
			print(paste('INFO: load previous GFIT result'))
			load(file=paste('Rdata/abc_result/gfit.',nSIM,'.rda',sep=''))
		}
		tGFIT<<-tGFIT
		for (Ne in Ne.vec){
			for (model in model.vec){
				pdf(paste('out_figures/Ne',Ne,'/gfit.',model,'.summary.pdf',sep=''))
				plot(tGFIT[[Ne]][[model]], main=paste(model,': p-value=',summary(tGFIT[[Ne]][[model]])$pvalue,sep=''))
				dev.off()
			}
		}
	} else print('WARNING: goodness of fit analyses disabled, please check bool.gfit')
}
## Summarize the result
in_summarize_accepted_distances=function(dist.list,Ne){
	df.dist=data.frame(do.call(rbind,dist.list))
	rownames(df.dist)=NULL
	df.dist$model=factor(df.dist$model,levels=model.vec)
	theme.common=theme(
		panel.background=element_rect(fill='white',color='black',size=0.5),
		text=element_text(size=7,hjust=.5,vjust=.5),
		line=element_line(size=0.5),
		panel.border=element_rect(color='black',linetype='solid',fill=NA),
		panel.grid.major=element_line(color='grey90',size=0.25),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=7,hjust=.5,vjust=.5),	
		plot.title=element_text(size=7,face='bold.italic'),
		plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		axis.line=element_blank(),
		axis.ticks=element_line(size=0.25),
		strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
		strip.background=element_rect(fill="grey90",color=NA,linetype='solid',size=0.5),
		legend.background=element_rect(fill='white',size=0.5),
		legend.margin=margin(unit(0.1,'cm')),
		legend.key=element_rect(color='black'),
		legend.key.size=unit(0.2,'lines'),
		legend.position='none')
		
	fig.dist=ggplot(df.dist)+
		geom_boxplot(aes(x=tolerance,y=accepted.dist,fill=model),width=0.8,color='black',alpha=0.8,lwd=0.8,outlier.shape=4,outlier.stroke=0.5)+
		scale_fill_manual(values=color.model)+
		labs(x='tolerance',y='accepted distance')+
		theme.common+
		facet_wrap(~model,ncol=4)
	ggsave(paste('out_figures/Ne',Ne,'/accepted_distance.pdf',sep=''),fig.dist,units='mm',height=80,width=150,dpi=600)
	out.list=list()
	i=1
	for (model in model.vec){
		for (tol in tol.vec){
			temp=df.dist[which(df.dist$model==model & df.dist$tolerance==tol),]
			out.list[[i]]=c(model=model,tolerance=tol,mean=round(mean(temp$accepted.dist),3),stdev=round(sd(temp$accepted.dist),3))
			i=i+1
		}
	}
	write.table(file=paste('out_tables/Ne',Ne,'/dist_summary.txt',sep=''),do.call(rbind,out.list),quote=F,col.names=T,row.names=F,sep='\t')
}
in_visualize_sumstat_plot=function(sumstat.list,Ne){
	theme.common=theme(
		panel.background=element_rect(fill='white',color='black',size=0.5),
		text=element_text(size=7,hjust=.5,vjust=.5),
		line=element_line(size=0.5),
		panel.border=element_rect(color='black',linetype='solid',fill=NA),
		panel.grid.major=element_line(color='grey90',size=0.25),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=7,hjust=.5,vjust=.5),	
		plot.title=element_text(size=7,face='bold.italic'),
		plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		axis.line=element_blank(),
		axis.ticks=element_line(size=0.25),
		strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
		strip.background=element_rect(fill="grey90",color=NA,linetype='solid',size=0.5),
		legend.background=element_rect(fill='white',size=0.5),
		legend.margin=margin(unit(0.1,'cm')),
		legend.key=element_rect(color='black'),
		legend.key.size=unit(0.2,'lines'),
		legend.position='none')
		
	for (tol in tol.vec){
		df.sumstat=sumstat.list[[tol]]
		rownames(df.sumstat)=NULL
		fig.sumstat=ggplot(df.sumstat)+
			geom_histogram(aes(x=value,y=..count..,fill=model),alpha=0.8,color=NA)+
			geom_vline(xintercept=obs.sumstat,color='black', size=0.5,linetype='longdash')+
			scale_fill_manual(values=color.model)+
			#stat_density_2d(geom='polygon',aes(alpha=..level..),color=NA,fill="dodgerblue2",size=0.5,bins=5)+
			#labs(x='GB',y='MN')+
			scale_x_continuous(expand=c(0,0.1))+
			scale_y_continuous(expand=c(0,0.1))+
			theme.common+
			facet_grid(model~variable)
		ggsave(paste('out_figures/Ne',Ne,'/sumstat_',tol,'.pdf',sep=''),fig.sumstat,units='mm',height=100,width=150,dpi=600)
	}
}
in_visualize_posterior_plot=function(posterior.list,Ne){
	for (tol in tol.vec){
		if (TRUE){
			temp=do.call(rbind,posterior.list)
			df.posterior=data.frame(type=rep(c('5% HPDI','mean','95% HPDI'),nrow(temp)/3),temp,stringsAsFactors=F)
			df.posterior$model=factor(df.posterior$model,levels=model.vec)
			rownames(df.posterior)=NULL
			df.posterior=data.frame(df.posterior[order(df.posterior$model),])
			for (i in 4:ncol(df.posterior)){
				df.posterior[,i]=round(as.numeric(df.posterior[,i]),3)
			}
			write.table(file=paste('out_tables/Ne',Ne,'/posterior_summary.txt',sep=''),df.posterior,quote=F,col.names=T,row.names=F,sep='\t')
			params.max=c(max_Fzm,max_Fzf,max_MF,max_Smn,max_Szs)
			params.min=c(0,0,min_MF,0,0)
			names(params.max)=colnames(df.posterior)[4:ncol(df.posterior)]
			names(params.min)=colnames(df.posterior)[4:ncol(df.posterior)]
		}
		if (TRUE) {
			fig.blank=ggplot()+geom_blank(aes(1,1))+
				theme(plot.background=element_blank(), 
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(), 
				panel.border=element_blank(),
				panel.background=element_blank(),
				plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"),
				axis.title=element_blank(),
				axis.text=element_blank(),
				axis.ticks=element_blank(),
				axis.line=element_blank(),
				legend.position='none')
			theme.1d=theme(text=element_text(size=7,hjust=.5,vjust=.5),
				line=element_line(size=0.5),
				panel.background=element_rect(fill='white', color='black'),
				panel.border=element_rect(color='black',linetype='solid',fill=NA),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_line(size=0.25,color='grey90'),
				strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
				axis.title.y=element_blank(),
				axis.ticks=element_line(size=0.25),
				plot.margin=unit(c(0.5,0.1,0.5,0.1),"cm"),
				legend.background=element_rect(fill='white',size=0.5),
				legend.margin=margin(unit(0.1,'cm')),
				legend.key=element_rect(color='black'),
				legend.key.size=unit(1,'lines'),
				legend.position='none')	
			theme.2d.mid=theme(text=element_text(size=7,hjust=.5,vjust=.5),
				line=element_line(size=0.5),
				panel.background=element_rect(fill='white', color='black'),
				panel.border=element_rect(color='black',linetype='solid',fill=NA),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_line(size=0.25,color='grey90'),
				strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
				axis.text=element_blank(),
				axis.title=element_blank(),
				axis.ticks=element_line(size=0.25),
				plot.margin=unit(c(0,0.1,0,0.1),"cm"),
				axis.line=element_blank(),
				legend.background=element_rect(fill='white',size=0.5),
				legend.margin=margin(unit(0.1,'cm')),
				legend.key=element_rect(color='black'),
				legend.key.size=unit(1,'lines'),
				legend.position='none')
			theme.2d.x=theme(text=element_text(size=7,hjust=.5,vjust=.5),
				line=element_line(size=0.5),
				panel.background=element_rect(fill='white', color='black'),
				panel.border=element_rect(color='black',linetype='solid',fill=NA),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_line(size=0.25,color='grey90'),
				strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
				axis.text.y=element_blank(),
				axis.title.y=element_blank(),
				axis.ticks=element_line(size=0.25),
				plot.margin=unit(c(0,0.1,0,0.1),"cm"),
				axis.line=element_blank(),
				legend.background=element_rect(fill='white',size=0.5),
				legend.margin=margin(unit(0.1,'cm')),
				legend.key=element_rect(color='black'),
				legend.key.size=unit(1,'lines'),
				legend.position='none')	
			theme.2d.y=theme(text=element_text(size=7,hjust=.5,vjust=.5),
				line=element_line(size=0.5),
				panel.background=element_rect(fill='white', color='black'),
				panel.border=element_rect(color='black',linetype='solid',fill=NA),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_line(size=0.25,color='grey90'),
				strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
				axis.text.x=element_blank(),
				axis.title.x=element_blank(),
				axis.ticks=element_line(size=0.25),
				plot.margin=unit(c(0,0.1,0,0.1),"cm"),
				axis.line=element_blank(),
				legend.background=element_rect(fill='white',size=0.5),
				legend.margin=margin(unit(0.1,'cm')),
				legend.key=element_rect(color='black'),
				legend.key.size=unit(1,'lines'),
				legend.position='none')		
		}
		bool.exclude=TRUE
		bool.1dplot=TRUE
		bool.2dplot=TRUE
		for (model in model.vec){
			print(paste('INFO: visualize posterior distribution, model=',model,', tolerance=',tol,sep=''))
			model.abc=tABC[[Ne]][[model]][[tol]]
			df.params=data.frame(tPARAMS[[Ne]][[model]],accepted=model.abc$region)
			df.stat=df.posterior[which(df.posterior$model==model & df.posterior$tolerance==tol),]
			n=ncol(df.params)-1
			if (bool.exclude){
				temp.stat=matrix(as.numeric(unlist(df.stat[,4:ncol(df.stat)])),nrow=3)
				exclude.vec=c(1:ncol(df.params))[which(apply(temp.stat,2,sum)==0)]
				if (length(exclude.vec) != 0){
					n=n-length(exclude.vec)
					df.params=df.params[,-c(exclude.vec)]
					df.stat=df.stat[,-c(exclude.vec+3)]
				}
			}
			if (bool.1dplot){
				tFIG.1D=vector('list',length=n)
				for (i in 1:n){
					name.i=colnames(df.params)[i]
					column.i=sym(name.i)
					mean.i=as.numeric(df.stat[which(df.stat$type=='mean'),name.i])
					lower.i=as.numeric(df.stat[which(df.stat$type=='5% HPDI'),name.i])
					higher.i=as.numeric(df.stat[which(df.stat$type=='95% HPDI'),name.i])	
					tFIG.1D[[i]]=ggplot(df.params)+
						geom_density(aes(x=!!column.i,group=accepted,fill=accepted),alpha=0.7,color=NA)+
						geom_vline(xintercept=mean.i,color='black',linetype='solid',alpha=0.8,size=0.5)+
						geom_vline(xintercept=c(lower.i,higher.i),color='black',linetype='longdash',alpha=0.8,size=0.5)+
						labs(x=name.i,y='')+
						scale_x_continuous(expand=c(0.05,0))+
						scale_y_continuous(expand=c(0.05,0))+
						scale_fill_manual(values=c('FALSE'='grey60','TRUE'='dodgerblue2'))+
						theme.1d
				}
				plot.1d=plot_grid(plotlist=tFIG.1D,ncol=n)
				save_plot(paste('out_figures/Ne',Ne,'/posterior_1d.',model,'.',tol,'.pdf',sep=''),plot.1d,units='mm',base_height=50,base_width=200,dpi=600)
			}
			if (bool.2dplot){
				tFIG.2D=vector('list',length=(n-1)^2)
				for (i in 1:(n-1)){
					for (j in 2:n){
						index=(n-1)*(j-1)-i+1
						if (index<=0) next
						temp.df=df.params[which(df.params$accepted),c(i,j)]
						name.i=colnames(df.params)[i]
						name.j=colnames(df.params)[j]
						column.i=sym(name.i)
						column.j=sym(name.j)
						mean.i=as.numeric(df.stat[which(df.stat$type=='mean'),name.i])
						mean.j=as.numeric(df.stat[which(df.stat$type=='mean'),name.j])
						if (i<j){
							tFIG.2D[[index]]=ggplot(temp.df,aes(x=!!column.i, y=!!column.j))+
								stat_density_2d(geom="polygon", aes(alpha=..level../sum(..level..)),fill='dodgerblue2')+
								geom_vline(xintercept=mean.i,color='black',linetype='solid',alpha=0.8,size=0.5)+
								geom_hline(yintercept=mean.j,color='black',linetype='solid',alpha=0.8,size=0.5)+
								scale_x_continuous(expand=c(0.05,0),limits=c(params.min[name.i],params.max[name.i]))+
								scale_y_continuous(expand=c(0.05,0),position='right',limits=c(params.min[name.j],params.max[name.j]))+
								labs(x=name.i,y=name.j,title='')+
								theme.2d.mid
								#{if(index%%(n-1)==0) theme.2d.y else if(index>(n-1)*(n-2)) theme.2d.x else theme.2d.mid}
						} else {
							tFIG.2D[[index]]=fig.blank
						}
					}
				}
				plot.2d=plot_grid(plotlist=tFIG.2D,ncol=n-1,nrow=n-1)
				save_plot(paste('out_figures/Ne',Ne,'/posterior_2d.',model,'.',tol,'.pdf',sep=''),plot.2d,units='mm',base_height=150,base_width=150,dpi=600)
			}
		}
	}
}
in_visualize_model_selection=function(Ne){
	theme.common=theme(text=element_text(size=7,hjust=.5,vjust=.5),
		line=element_line(size=0.5),
		panel.background=element_rect(fill='white', color='black'),
		panel.border=element_rect(color='black',linetype='solid',fill=NA),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_line(size=0.25,color='grey90'),
		strip.text=element_text(color="black",size=10,margin=margin(1,1,1,1,'mm'),face="bold"),
		axis.ticks=element_line(size=0.25),
		plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
		legend.background=element_rect(fill='white',size=0.5),
		legend.margin=margin(unit(0.1,'cm')),
		legend.key=element_rect(color='black'),
		legend.key.size=unit(1,'lines'),
		legend.position='right')
	cv.list=list()
	recall.list=list()
	title.string=''
	for (tol in tol.vec){
		cv.modsel=tMODSEL[[Ne]][[tol]][['cvmodsel']]
		conf.mat=summary(cv.modsel)$conf.matrix[[1]]
		write.table(file=paste('out_tables/Ne',Ne,'/confusion_matrix.',tol,'.txt',sep=''),conf.mat,quote=F,col.names=T,row.names=F,sep='\t')
		temp.list=list()
		for (model in model.vec){
			row=conf.mat[model,]
			temp.list[[model]]=cbind(true.model=rep(model,length(row)),predict.model=names(row),frequency=row/sum(row))
		}
		cv.list[[tol]]=data.frame(tolerance=rep(tol,length(row)*length(model.vec)),do.call(rbind,temp.list))
		recall.list[[tol]]=round(sum(diag(conf.mat))/sum(conf.mat),3)
		modsel=tMODSEL[[Ne]][[tol]][['modsel']]
		bayes.factor=summary(modsel)$BayesF
		write.table(file=paste('out_tables/Ne',Ne,'/modelselection_',tol,'.txt',sep=''),round(bayes.factor,3),quote=F,col.names=T,row.names=T,sep='\t')
		title.string=paste(title.string,'\nrecall_',tol,'=',recall.list[[tol]],sep='')
	}
	df.conf=do.call(rbind,cv.list)
	df.conf$frequency=as.numeric(as.character(df.conf$frequency))
	rownames(df.conf)=NULL
	df.conf$true.model=factor(df.conf$true.model,levels=model.vec)
	df.conf$predict.model=factor(df.conf$predict.model,levels=model.vec)
	fig.conf=ggplot(df.conf,aes(x=predict.model,y=true.model))+
		geom_tile(aes(fill=frequency),color=NA,size=0.5)+
		geom_text(aes(label=frequency),color='black',size=1.5,hjust=.5,vjust=.5)+
		scale_fill_gradient(high="dodgerblue2", low="white", space="Lab",  na.value="black", guide="colourbar",name='')+
		scale_y_discrete(expand=c(0.05,0))+
		scale_x_discrete(expand=c(0.05,0))+
		labs(x='true model',y='predicted model',title=title.string)+
		coord_fixed(ratio=1)+
		theme.common+
		facet_wrap(~tolerance)
	ggsave(paste('out_figures/Ne',Ne,'/confusion_matrix.pdf',sep=''),fig.conf,units='mm',height=70,width=150,dpi=600)
}
summarize=function(){
	for (Ne in Ne.vec){
		dir.create(file.path(paste('out_figures/Ne',Ne,sep='')), showWarnings=FALSE)
		dir.create(file.path(paste('out_tables/Ne',Ne,sep='')), showWarnings=FALSE)
		dist.list=list()
		sumstat.list=list()
		posterior.list=list()
		for (tol in tol.vec){
			temp.dist=list()
			temp.sumstat=list()
			temp.posterior=list()
			for (model in model.vec){	
				model.abc=tABC[[Ne]][[model]][[tol]]
				dist.vec=model.abc$dist[model.abc$region]
				temp.dist[[model]]=data.frame(model=rep(model,length(dist.vec)),tolerance=rep(tol,length(dist.vec)),accepted.dist=as.numeric(dist.vec))
				temp.posterior[[model]]=cbind(model=rep(model,3),tolerance=rep(tol,3),summary(model.abc,intvl=0.9)[c(2,4,6),])
				temp.mat=tMATRIX[[Ne]][[model]][which(model.abc$region),]
				temp.melt=list()
				i=1
				for (column in 7:ncol(temp.mat)){
					temp.melt[[i]]=data.frame(model=rep(model,nrow(temp.mat)),tolerance=rep(tol,nrow(temp.mat)),variable=rep(colnames(temp.mat)[column],nrow(temp.mat)),value=as.numeric(as.character(temp.mat[,column])))
					i=i+1
				}
				temp.sumstat[[model]]=do.call(rbind,temp.melt)
			}
			sumstat.list[[tol]]=do.call(rbind,temp.sumstat)
			dist.list[[tol]]=do.call(rbind,temp.dist)
			posterior.list[[tol]]=do.call(rbind,temp.posterior)
		}
		## Distance
		in_summarize_accepted_distances(dist.list,Ne)
		## Sumstat plot
		in_visualize_sumstat_plot(sumstat.list,Ne)
		## Posterior inference
		in_visualize_posterior_plot(posterior.list,Ne)
		## Model selection
		in_visualize_model_selection(Ne)
	}
}

load_data()
run_abc()
run_model_selection()
summarize()
run_goodness_of_fit()
head(warnings(),20)
