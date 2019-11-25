####################################################################################
############################MADE BY TAEHYUNG KWON Nov 11 2019#######################
####################################################################################

##PRE-REQUISITE: python3, numpy
import sys, os, time
import numpy as np
import threading

##GLOBAL PARAMETERS
nSIM = 20000000
nTHREAD = 24
nSET = 100000
##SIMULATION PARAMETERS
nREP = 1
nM = 200	#Male popsize
nF = 800	#Female popsize
nGEN = 150	#Number of generation

##ABC PARAMETERS
GB_obs = 0.755080218
MT_obs = 0
TOL = 0.001	#Tolerence of ABC

####################################################################################
def run_simulation(nSET, THR_NUM):
	sema.acquire()
	vars = map(str,[THR_NUM,nSET,nREP,nM,nF,nGEN])
	start = time.time()
	print('RUN STARTED: {0} POPULATION, {1}/{2}'.format(nSET,THR_NUM,nSIM/nSET))
	os.system('Rscript run_sim.R '+' '.join(vars))
	end = time.time()
	os.system('echo \'SET{0}\t{1}\t{2}\t{3}sec\' >> ./time_log_nSIM{4}.txt'.format(THR_NUM,nSET,nTHREAD,round(end-start,0),nSIM))
	sema.release()

def simulate():
	ths = []
	global sema
	sema = threading.Semaphore(value=nTHREAD)
	i=1
	os.system('mkdir Rdata')
	os.system('echo \'SET\tnSET\tnTHREAD\tTIME\' > ./time_log_nSIM{0}.txt'.format(nSIM))#log file
	while (i <= nSIM/nSET):
		th = threading.Thread(target=run_simulation, args=(nSET,i,))
		th.daemon = True
		ths.append(th)
		th.start()
		i += 1

	for thread in ths:
		thread.join()
####################################################################################
def summarize():
	import glob
	os.system('mkdir results')
	total = len(glob.glob('Rdata/set*.rda'))
	os.system('Rscript run_abc.R {0} {1} {2} {3} {4}'.format(total,nSIM,GB_obs,MT_obs,TOL))
		
####################################################################################
##EXECUTE
#simulate()
summarize()
