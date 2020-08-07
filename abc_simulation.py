####################################################################################
############################MADE BY TAEHYUNG KWON Nov 11 2019#######################
########################CORRECTED BY TAEHYUNG KWON Feb 18 2020######################
####################################################################################

##PRE-REQUISITE: python3, numpy
import sys, os, time
import numpy as np
import threading

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

####################################################################################
def run_simulation(nSET,THR_NUM,timestr):
	sema.acquire()
	vars = map(str,[THR_NUM,nSET,nGEN,nPOP,max_Fzm,max_Fzf,max_MF,max_Smn,max_Szs])
	
	start = time.time()
	print('RUN STARTED: {0} REPLICATES, {1}/{2}'.format(nSET,THR_NUM,nSIM/nSET))
	os.system('Rscript run_sim.R '+' '.join(vars))
	
	end = time.time()
	os.system('echo \'SET{0}\t{1}\t{2}\t{3}sec\' >> ./time_{4}.log.txt'.format(THR_NUM,nSET,nTHREAD,round(end-start,0),timestr))	
	sema.release()
	
def simulate():
	ths = []
	global sema
	sema = threading.Semaphore(value=nTHREAD)
	i=1
	os.system('mkdir -p Rdata')
	timestr=time.strftime('%m-%d-%H:%M')
	os.system('echo \'SET\tnSET\tnTHREAD\tTIME\' > ./time_{0}.log.txt'.format(timestr))
	while (i <= nSIM/nSET):
		th = threading.Thread(target=run_simulation, args=(nSET,i,timestr,))
		th.daemon = True
		ths.append(th)
		th.start()
		i += 1

	for thread in ths:
		thread.join()

####################################################################################
def summarize():
	import glob
	os.system('mkdir -p results figures')
	total = len(glob.glob('Rdata/a.set*.rda'))
	os.system('Rscript run_abc.R {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(total,nSIM,GB_obs,MT_obs,Y_obs,TOL,max_Fzm,max_Fzf,max_MF,max_Smn,max_Szs))
		
####################################################################################
def main():
	if (len(sys.argv)==1):
		print('ERROR:\tPLEASE INPUT OPTION\n\t--sim\t\tGENERATE SIMULATES FOR ABC\n\t--sum\t\tPERFORM ABC AND VISUALIZE THE RESULT')
	elif (sys.argv[1] == '--sim'):
		print('INFO:\tSIMULATION STARTED\n{0}'.format('\t\n'.join(['nSIM='+str(nSIM),'nTHREAD='+str(nTHREAD),'nSET='+str(nSET),'Population size='+str(nPOP),'Generation number='+str(nGEN)])))
		simulate()
	elif (sys.argv[1] == '--sum'):
		print('ABC STARTED:\n{0}'.format('\tTolerance='+str(TOL)))
		summarize()
	else:
		print('ERROR:\tWRONG OPTION SPECIFIED. PLEASE TRY "--sim" OR "--sum"')

main()
