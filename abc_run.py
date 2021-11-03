### MADE BY TAEHYUNG KWON
## Requisites: python3, R
import sys, os, time, glob, argparse, multiprocessing
from argparse import RawTextHelpFormatter

## GLOBAL PARAMETERS
nSIM=10000000
nREP=10000

## SIMULATION PARAMETERS
Ne_1=500
Ne_2=2000
Ne_3=5000
ngen=110
max_Fzm=1
max_Fzf=0.5
min_MF=0.05
max_MF=0.5
max_Smn=0.2
max_Szs=100

## ABC PARAMETERS
GB_mean=0.75508
MT_mean=0
Y_mean=1

def __run__simulation(nREP,PRC,timestr):
	sema.acquire()
	tVARS=' '.join(list(map(str,[PRC,nREP])))
	start=time.time()
	print(f'INFO: run replicate ({PRC}/{int(nSIM/nREP)})')
	os.system(f'Rscript {workdir}/abc_simulate.R {tVARS}')
	end=time.time()
	os.system('echo \'{0}\t{1}\t{2}\t{3}sec\' >> {4}/time_{5}.log.txt'.format(PRC,nREP,nCPU,round(end-start,0),workdir,timestr))
	sema.release()
	
def simulate():
	tPROCESS=[]
	global sema
	sema=multiprocessing.Semaphore(value=nCPU)
	PRC=1
	timestr=time.strftime('%m-%d-%H:%M')
	os.system(f'echo \'replicate_number\treplicate_size\ttotal_CPU\ttime\' > {workdir}/time_{timestr}.log.txt')
	print(f'INFO:\tsimulation ready\n\ttotal simulation replicates={nSIM}\n\treplicate size={nREP}')
	while (PRC <= nSIM/nREP):
		if len(glob.glob(f'{workdir}/Rdata/*.Ne{Ne_1}.rep{PRC}.rda'))==4 and len(glob.glob(f'{workdir}/Rdata/*.Ne{Ne_2}.rep{PRC}.rda'))==4 and len(glob.glob(f'{workdir}/Rdata/*.Ne{Ne_3}.rep{PRC}.rda'))==4:
			print(f'WARNING: skip replicate ({PRC}/{int(nSIM/nREP)})')
		else:
			process=multiprocessing.Process(target=__run__simulation, args=(nREP,PRC,timestr,))
			process.daemon=True
			tPROCESS.append(process)
			process.start()
		PRC += 1
	for process in tPROCESS:
		process.join()

def summarize():
	os.system(f'mkdir -p {workdir}/out_tables {workdir}/out_figures')
	total=len(glob.glob(f'{workdir}/Rdata/bothsel.Ne{Ne_3}.rep*.rda'))
	os.system(f'Rscript {workdir}/abc_summarize.R {total} {nREP}')

def save_parameter():
	os.system(f'mkdir -p {workdir}/Rdata')
	fparam_sim=f'{workdir}/Rdata/params.sim.txt'
	os.system('echo "Ne_1\t{}\nNe_2\t{}\nNe_3\t{}\nngen\t{}\nmax_Fzm\t{}\nmax_Fzf\t{}\nmin_MF\t{}\nmax_MF\t{}\nmax_Smn\t{}\nmax_Szs\t{}" > {}'.format(
		Ne_1,
		Ne_2,
		Ne_3,
		ngen,
		max_Fzm,
		max_Fzf,
		min_MF,
		max_MF,
		max_Smn,
		max_Szs,
		fparam_sim
		))
	fparam_abc=f'{workdir}/Rdata/params.abc.txt'
	os.system('echo "GB_mean\t{}\nMT_mean\t{}\nY_mean\t{}" > {}'.format(
		GB_mean,
		MT_mean,
		Y_mean,
		fparam_abc
		))

def __main__():
	parser=argparse.ArgumentParser(description='ABC for demographic simulation by Taehyung Kwon', formatter_class=RawTextHelpFormatter)
	parser.add_argument('function', type=str,
		choices=['sim','sum'],
		help='perform forward simulation: "sim"\n\
			summarize simulation result: "sum"')
	parser.add_argument('-p', '--cpu', type=int, default=1, help='number of processors')
	parser.add_argument('-d', '--dir', type=str, default=os.getcwd(), help='working directory')
	args=parser.parse_args()
	global nCPU, workdir
	nCPU=args.cpu
	workdir=args.dir
	if workdir.endswith('/'):
		workdir=workdir[:-1]
	save_parameter()
	if args.function == 'sim':
		simulate()
	elif args.function == 'sum':
		summarize()
		
__main__()