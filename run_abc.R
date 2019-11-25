rm(list=ls())
####################################################################################
############################MADE BY TAEHYUNG KWON Nov 11 2019#######################
####################################################################################

##REQUISITE FOR ABC: abc
##REQUISITE FOR PLOT: ggplot2, plyr, cowplot, egg, ggExtra

abc_plot=function(nSIM,PARAMS,abc.result,max.FM0,max.PZ,max.MNS,max.ZMS,TOL){
	library(ggplot2)
	library(plyr)
	library(cowplot)
	library(egg)
	library(ggExtra)
	print('ABC plots generating ...')
	
	## ABC SUMMARY TABLE
	abc_med=summary(abc.result,print=FALSE)[3,]				#median
	abc_q05=summary(abc.result,print=FALSE,intvl=0.9)[2,]	#5% quantile
	abc_q95=summary(abc.result,print=FALSE,intvl=0.9)[6,]	#95% quantile
	
	sum.stat=rbind(c('mean.GB','sd.GB','mean.MT','sd.MT','mean.MN','sd.MN'),
	cbind(mean(abc.result$ss[,1]),sd(abc.result$ss[,1]),mean(abc.result$ss[,2]),sd(abc.result$ss[,2]),mean(mDATA[which(abc.result$region),7]),sd(mDATA[which(abc.result$region),7])))
		
	write.table(file=paste('results/params_summary.tol',TOL,'.txt',sep=''),rbind(abc_med,abc_q05,abc_q95),quote=F,col.names=T,row.names=T,sep='\t')
	write.table(file=paste('results/sumstat_summary.tol',TOL,'.txt',sep=''),sum.stat,quote=F,col.names=F,row.names=F,sep='\t')
	
	## HEATMAP FOR ACCEPTANCE
	max.params=c(max.FM0,max.PZ,max.MNS,max.ZMS)
	break.params=c(max.FM0/10,max.PZ/10,max.MNS/10,max.ZMS/10)
	ALL=data.frame(cut(PARAMS[,1],breaks=seq(0,max.params[1],break.params[1]),right=FALSE),
		cut(PARAMS[,2],breaks=seq(0,max.params[2],break.params[2]),right=FALSE),
		cut(PARAMS[,3],breaks=seq(0,max.params[3],break.params[3]),right=FALSE),
		cut(PARAMS[,4],breaks=seq(0,max.params[4],break.params[4]),right=FALSE),
		abc.result$region)
	colnames(ALL)=c('FM0','PZ','MNS','ZMS','ACC')
	ALL[which(ALL$ACC==TRUE),5]=1
	
	temp=ddply(ALL,c('FM0','PZ','MNS','ZMS'),numcolwise(sum))
	counter=ddply(ALL,c('FM0','PZ','MNS','ZMS'),nrow)
	temp$ACC=round(temp$ACC/counter$V1,3)
	temp$ZMS=factor(temp$ZMS, levels=rev(levels(temp$ZMS)))
	fig.heatmap=ggplot(temp,aes(as.factor(FM0),as.factor(PZ)))+
		geom_tile(aes(fill=ACC),color=NA,size=0.1)+
		scale_fill_gradient(name='Acceptance rates',limits=c(0,max(temp$ACC)),low="black", high='white')+
		scale_y_discrete(expand=c(0,0),breaks=c(0,1),labels=c('0','1'))+
		scale_x_discrete(expand=c(0,0),breaks=c(0,1),labels=c('0','1'))+
		labs(x='',y='',title='')+
		coord_fixed(ratio=1)+
		theme(panel.background=element_rect(fill='white', colour='white',size=0.1),
			text=element_text(size=6,hjust=.5,vjust=.5),
			line=element_line(size=0.25),
			panel.border=element_rect(color='white',linetype='solid',size=0.1,fill=NA),
			panel.grid.minor=element_blank(),
			panel.grid.major=element_blank(),
			panel.spacing=unit(0.10,'lines'),
			axis.text.y=element_text(size=8,hjust=.5,vjust=.5),
			axis.text.x=element_text(size=8,hjust=.5,vjust=.5,angle=90),
			plot.title=element_text(size=8,face='bold.italic'),
			plot.margin=unit(c(1,0,1,0.5),"lines"),
			axis.line=element_blank(),
			strip.text=element_text(colour="white",size=6,margin=margin(0.5,0.1,0.5,0.1,'mm'),face="bold"),
			strip.background=element_rect(fill="#2E3842",color="white",linetype='solid',size=0.1),
			legend.background=element_rect(fill='white',size=0.2),
			legend.margin=margin(unit(0.1,'cm')),
			legend.key=element_rect(color='black'),
			legend.key.size=unit(1,'lines'),
			legend.position='bottom')+
		facet_grid(ZMS~MNS)
	
	ggsave(paste('results/heatmap_tol',TOL,'.png',sep=''),fig.heatmap,units='mm',height=150,width=150,dpi=600)
	ggsave(paste('results/heatmap_tol',TOL,'.pdf',sep=''),fig.heatmap,units='mm',height=150,width=150,dpi=300)	
	
	
	## POSTERIOR PROB
	break.params=c(max.FM0/50,max.PZ/50,max.MNS/50,max.ZMS/50)
	
	FORHIST=data.frame(cut(PARAMS[,1],labels=seq(break.params[1]*0.5,max.params[1]-break.params[1]*0.5,break.params[1]),breaks=seq(0,max.params[1],break.params[1]),right=FALSE),
		cut(PARAMS[,2],labels=seq(break.params[2]*0.5,max.params[2]-break.params[2]*0.5,break.params[2]),breaks=seq(0,max.params[2],break.params[2]),right=FALSE),
		cut(PARAMS[,3],labels=seq(break.params[3]*0.5,max.params[3],break.params[3]),breaks=seq(0,max.params[3],break.params[3]),right=FALSE),
		cut(PARAMS[,4],labels=seq(break.params[4]*0.5,max.params[4],break.params[4]),breaks=seq(0,max.params[4],break.params[4]),right=FALSE),
		abc.result$region)
	
	colnames(FORHIST)=c('FM0','PZ','MNS','ZMS','ACC')
	
	FORHIST[which(FORHIST$ACC==TRUE),5]=1
	fig_one=list()
	for (i in c(1,2,3,4)){
		
		name.i=colnames(FORHIST)[i]
		temp.acc=ddply(FORHIST[,c(i,5)],c(name.i),numcolwise(sum))
		temp.pri=ddply(FORHIST[,c(i,5)],c(name.i),nrow)
		temp.hist=data.frame(var1=as.double(as.character(temp.acc[,1])),ACC=as.double(temp.acc$ACC/(nSIM*TOL)),PRIOR=as.double(temp.pri$V1/nSIM))
		
		fig.pp=ggplot(temp.hist)+
			geom_bar(aes(x=as.double(var1), y=PRIOR),stat='identity',fill='black',alpha=0.2,position=position_dodge(1))+
			geom_bar(aes(x=as.double(var1), y=ACC),fill='red',stat='identity',alpha=0.5,position=position_dodge(1))+
			geom_vline(xintercept=abc_med[i],color='black',linetype='solid')+
			geom_vline(xintercept=c(abc_q05[i],abc_q95[i]),color='black',linetype='longdash')+
			scale_x_continuous(expand=c(0,0.001))+
			scale_y_continuous(expand=c(0,0.001),limits=c(0,0.07))+
			labs(x='',y='',title='')+
			theme(panel.background=element_rect(fill='white', colour='black',size=0.2),
				text=element_text(size=6,hjust=.5,vjust=.5),
				line=element_line(size=0.25),
				panel.border=element_rect(color='black',linetype='solid',size=0.2,fill=NA),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_line(size=0.05,color='grey'),
				panel.spacing=unit(0.15,'lines'),
				axis.text.y=element_text(size=8,hjust=.5,vjust=.5),
				axis.text.x=element_text(size=8,hjust=.5,vjust=.5,angle=90),	
				plot.title=element_text(size=8,face='bold.italic'),
				plot.margin=unit(c(0,0.5,0,0.5),"lines"),
				axis.line=element_blank(),
				strip.text=element_text(colour="white",size=8,margin=margin(0,0,0,0,'mm'),face="bold"),
				strip.background=element_rect(fill="#1E403C",color="black",linetype='solid',size=0.2),
				legend.background=element_rect(fill='white',size=0.2),
				legend.margin=margin(unit(0.1,'cm')),
				legend.key=element_rect(color='black'),
				legend.key.size=unit(1,'lines'),
				legend.position='bottom')

		fig_one[[i]]=fig.pp
		}
	finalplot_one=ggarrange(fig_one[[1]],fig_one[[2]],fig_one[[3]],fig_one[[4]], nrow=4)
	save_plot(paste('results/ppdist_tol',TOL,'.pdf',sep=''),finalplot_one,units='mm',base_height=150,base_width=50,dpi=600)
	
	
	## SUMSTAT PLOT
	library(ggExtra)
	
	df.ss=data.frame(abc.result$ss)
	fig.ss=ggplot(df.ss)+
		geom_point(aes(x=df.ss$GB.150,y=MT.150),color='black',size=0.01)+
		geom_point(x=GB_obs,y=MT_obs,color='red', size=5, shape=18)+
		labs(x='GB',y='MT')+
		theme(panel.background=element_rect(fill='white', colour='black',size=0.2),
			text=element_text(size=6,hjust=.5,vjust=.5),
			line=element_line(size=0.25),
			panel.border=element_rect(color='black',linetype='solid',fill=NA),
			panel.grid.major=element_line(color='grey'),
			panel.grid.minor=element_line(color='grey'),
			panel.spacing=unit(0.15,'lines'),
			axis.text.y=element_text(size=8,hjust=.5,vjust=.5),
			axis.text.x=element_text(size=8,hjust=.5,vjust=.5,angle=90),	
			plot.title=element_text(size=8,face='bold.italic'),
			plot.margin=unit(c(0.5,0.5,0.5,0.5),"lines"),
			axis.line=element_blank(),
			strip.text=element_text(colour="white",size=8,margin=margin(0,0,0,0,'mm'),face="bold"),
			strip.background=element_rect(fill="#1E403C",color="black",linetype='solid',size=0.2),
			legend.background=element_rect(fill='white',size=0.2),
			legend.margin=margin(unit(0.1,'cm')),
			legend.key=element_rect(color='black'),
			legend.key.size=unit(1,'lines'),
			legend.position='bottom')
	
	fig.ss.marg = ggMarginal(fig.ss, type="histogram",margins='both',size=5, color=NA,fill='black')
	ggsave(paste('results/sumstat_tol',TOL,'.png',sep=''),fig.ss.marg,units='mm',height=170,width=170,dpi=600)
		
	warnings()
	print('ABC plots generated')
}

execute_abc=function(){
	library(abc)

	## TEMP PARAMETERS
	#GB_obs=0.755080218
	#MT_obs=0
	#TOTAL.SET=500
	#nSIM=20000000
	#TOL=0.01

	## LIMITS OF PARAMETERS
	max.FM0=1
	max.PZ=1
	max.MNS=0.2
	max.ZMS=50
	nVAL=100
	set=1
	tDATA=list()
	print('simulation set loading ...')
	while (set <= TOTAL.SET){ 
		infile=paste('Rdata/set',set,'.rda',sep='')
		if (file.exists(infile)){
			load(file=infile)
			if (ncol(tSUMMARY)>50){
				tDATA[[set]]=tSUMMARY[,c(1:4,53,54,55)]
			} else {
				tDATA[[set]]=tSUMMARY
			}
		}
		set=set + 1
	}
	print('all simulation set loaded')
	mDATA=do.call('rbind',tDATA)
	OBS.SS=c(GB_obs,MT_obs)
	PARAMS=mDATA[1:nSIM,1:4]
	SIM.SS=mDATA[1:nSIM,c(5,6)]
	
	## ABC - takes few minutes
	TOL=0.01
	if (file.exists(paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))==FALSE){
		start.time=Sys.time()
		print('ABC performing ...')
		abc.result.01=abc(target=OBS.SS, param=PARAMS, sumstat=SIM.SS,tol=TOL,method="rejection")
		end.time=Sys.time()
		save(abc.result, file=paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))
		print(paste('ABC performed: ',end.time-start.time,sep=''))
	} else {
		load(file=paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))
		print(paste('ABC skipped: ABC result exists!'))
	}

	TOL=0.001
	if (file.exists(paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))==FALSE){
		start.time=Sys.time()
		print('ABC performing ...')
		abc.result.001=abc(target=OBS.SS, param=PARAMS, sumstat=SIM.SS,tol=TOL,method="rejection")
		end.time=Sys.time()
		save(abc.result, file=paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))
		print(paste('ABC performed: ',end.time-start.time,sep=''))
	} else {
		load(file=paste('Rdata/ABC_RES.',nSIM,'_',TOL,'.rda',sep=''))
		print(paste('ABC skipped: ABC result exists!'))
	}
	
	summary(abc.result.01$ss)
	summary(abc.result.001$ss)
	
	## GENERATE ABC PLOT
	quit()
	abc_plot(nSIM,PARAMS,abc.result.01,max.FM0,max.PZ,max.MNS,max.ZMS,0.01)
	abc_plot(nSIM,PARAMS,abc.result.001,max.FM0,max.PZ,max.MNS,max.ZMS,0.001)
	
}	
####################################################################################
args=commandArgs(trailingOnly=TRUE)
TOTAL.SET <<- as.integer(args[1])
nSIM <<- as.integer(args[2])
GB_obs <<- as.double(args[3])
MT_obs <<- as.double(args[4])
TOL <<- as.double(args[5])

execute_abc()
