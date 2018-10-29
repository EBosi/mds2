#!/usr/bin/python

"""
MeDuSa: A draft genome scaffolder that uses multiple 
reference genomes in a graph-based approach
"""

#############
# Imports   #
#############
import netcon_mummer
from medusa_lib import Scaffolder
from medusa_lib import *
#from graph_lib import *    
import logging,sys,os
from optparse import OptionParser,OptionGroup
import tempfile
import time
import multiprocessing_lib

#################
# Opt parsing   #
#################
if __name__ == "__main__":

	logo = """

 /$$      /$$               /$$                            /$$$$$$ 
| $$$    /$$$              | $$                           /$$__  $$
| $$$$  /$$$$ /$$$$$$  /$$$$$$$/$$   /$$ /$$$$$$$ /$$$$$$|__/  \ $$
| $$ $$/$$ $$/$$__  $$/$$__  $| $$  | $$/$$_____/|____  $$ /$$$$$$/
| $$  $$$| $| $$$$$$$| $$  | $| $$  | $|  $$$$$$  /$$$$$$$/$$____/ 
| $$\  $ | $| $$_____| $$  | $| $$  | $$\____  $$/$$__  $| $$      
| $$ \/  | $|  $$$$$$|  $$$$$$|  $$$$$$//$$$$$$$|  $$$$$$| $$$$$$$$
|__/     |__/\_______/\_______/\______/|_______/ \_______|________/
                                                                   

"""

	print logo

	usage=""" 



python %prog [options]
	"""

	def outGraphCallback1(option, opt_str, value, parser):
		if value == None: value = ''
		else: value += '_'
		setattr(parser.values, option.dest, value)	

	parser = OptionParser(usage=usage)
	
	# mandatory
	group1 = OptionGroup(parser, "Mandatory Arguments")
	group1.add_option("-i", "--input", dest="target",
	                  help=" target genome to be scaffolded ", metavar="FILE")
	group1.add_option("-f", "--files", dest="comparison_dir",
	                  help=" DIR where the comparison genomes are located ", metavar="DIR")
	parser.add_option_group(group1)
	
	# optional
	group2 = OptionGroup(parser, "Optional Arguments")
	group2.add_option("-o", "--output_dir", dest="output_dir",default=None,
	                  help="""
write output to DIR. If unset Medusa will generate a new directory formatted using local time as it follows: ./DDMMYYYY.HHMMSS/""", metavar="DIR")
	group2.add_option("--verboseMUMmer", dest="verbose", action="store_true",default=False,
	                  help=" print to STDOUT the information given by MUMmer ")
	group2.add_option("-r", "--random", dest="random", action="store_true",default=False,
	                  help="allow for random choice when best edges have the same score. This might lead to variable results!")
	group2.add_option("-w", "--weight", dest="weightScheme2", action="store_true",default=False,
	                  help=" allow for sequence similarity based weighting scheme. May lead to better results.")
	group2.add_option("-s", "--skipMapping", dest="skipMap",default=None,
	                  help=" Skip the mapping step. MUMmer output files will be searched in DIR. Option -m will be ignored ",metavar = "DIR")
	group2.add_option("-g", "--graph", dest="outGraph", metavar="PREFIX",action="callback",callback=outGraphCallback1,
	                  help=""" write scaffolding network and path covers using the prefix FILE in the output directory. 
Without any prefix, names will have a generic one """)
	group2.add_option("-d", "--distance", dest="distance", action="store_true",default=False,
	                  help=""" allow for the estimation of distance between pairs of contigs 
								based ond the reference genome(s): in this case the scaffolded contigs
								will be separated by a number of N characters equal to the estimate. 
								The estimated distances are also saved in the <targetGenome>_distanceTable file.
								By default, the scaffolded contigs are separated by 100 Ns. """)
	group2.add_option("-c","--cleanUp",dest="cleanUp",action="store_true",default=False, 
						help="clean temporary files (mummer files and graphs)")
	group2.add_option("-t","--threads",dest="threads",default=1,help="number of threads for multiprocessing")
	group2.add_option("-v","--verbosity",dest="logLevel",type=int,default=2,help="set the logger verbosity level to 0,1,2 or 3 (no log,warning,info,debug)")
	group2.add_option("--inputGraph",dest="inputGraph",default=None,help="""
resume analysis from a previous scaffolding graph. The option `-f` will be ignored.""")	


	parser.add_option_group(group2)
	
	(options, args) = parser.parse_args()
	if not options.target or not options.comparison_dir:
	    parser.print_help()
	    parser.error('Mandatory Arguments missing')

	
	##########
	# Logger #
	##########
	
	logging_level = [logging.ERROR,logging.WARNING,logging.INFO,logging.DEBUG][options.logLevel]
	logger = logging.getLogger(__name__)
	logging.basicConfig(level=logging_level,format='[%(asctime)s] %(levelname)-8s %(message)s',datefmt='%m-%d %H:%M:%S')
	
	
	
	########
	# Main #
	########
	wd = options.output_dir
	if wd == None: wd = time.strftime("./%d%m%y.%H%M%S/")
	target,comparison_dir,wd = map(os.path.abspath,[options.target,options.comparison_dir,wd])
	
	
	logging.info(greenText('Phase 0: checking the inputs...\n'))
	
	logging.info('Working directory is [%s]...' %wd)
	checkExistence(target);checkExistence(comparison_dir)
	logging.info(greenText('Input files found'))
	coords = []
	
	coordsDir = getMummerOutDir(wd)
	
	
	## do MUMmer for each pair
	logging.info(greenText('Phase 1: mapping contigs to reference genome(s)\n'))
	if options.skipMap == None and options.inputGraph == None:
	
		comparisons = [os.path.join(comparison_dir,f) for f in os.listdir(comparison_dir) if os.path.join(comparison_dir,f) != target]
		n = len(comparisons)
		if n > 1: logging.info('A total of %s reference genomes have been found! Mapping the contigs...\n' %str(n))
		else:
			logging.error(redText('No files have been found in %s... please change your input and retry' %comparison_dir))
			sys.exit()
		os.chdir(coordsDir)
		iterable_mummer_args = [(
						(target,c),
						{'verbose':options.verbose,'outputDir':options.output_dir,
						'threads':options.threads,'logging_info':'File number %s: %s' %(i+1,c),
						'fxn':runMummer}
					) for i,c in enumerate(comparisons)]

		coords = multiprocessing_lib.poolWrapper(int(options.threads),iterable_mummer_args)
					
	else:
		if options.inputGraph == None:
			coordsDir = os.path.abspath(options.skipMap)
			logging.info(greenText('Option -s (--skipMap) selected. MUMmer .coords files will be searched in %s' %coordsDir))
			checkExistence(coordsDir)
			coords = [os.path.join(coordsDir,f) for f in os.listdir(coordsDir) if f.endswith('.coords')]
		else:
			if options.skipMap != None and options.inputGraph != None:
				logging.info(yellowText('\nOptions -s (--skipMap) and --inputGraph were selected... --skipMap will be ignored!'))
			logging.info(greenText('\nOption --inputGraph selected... will directly analyse the scaffolding graph %s' %options.inputGraph))
			providedGraph = os.path.abspath(options.inputGraph)
			checkExistence(providedGraph)

	## from MUMmer to network
	readNet = True
	try: Scaffolding_graph = providedGraph
	except:
		readNet = False
		print
		logging.info(greenText('Phase 2: parse mummer output(s) and create a scaffolding graph'))
		logging.warning(yellowText('This step is very time-consuming...'))
	
		Scaffolding_graph = netcon_mummer.initialize_graph(target)
		if len(coords) > 1: logging.info('A total of %s .coords files have been found! Adding network edges...\n' %str(len(coords)))
		else:
			logging.error(redText('No files have been found in %s... please change your input and retry' %comparison_dir))
			sys.exit()
		netcon_mummer.coords2graph(coords,Scaffolding_graph,distanceEstimation=options.distance,altWeightScheme=options.weightScheme2)
	# from network 2 scaffold
	logging.info(greenText('Phase 3: from network to scaffolds\n'))
	scaffolder = Scaffolder()
	outScaffoldsFileName = os.path.abspath(os.path.join(wd,'scaffolds'))
	scaffolder.job(Scaffolding_graph,outScaffoldsFileName,target,readNet=readNet,threads=options.threads,distanceEstimation=options.distance)
	
	# export results
	if options.cleanUp:
		logging.info(greenText('Optional: cleaning up output directory'))
		cleanUp(coordsDir)
	if options.outGraph != None:
		logging.info(greenText('Optional: producing output graphs'))
		prefix = os.path.join(wd,options.outGraph)
		scaffolder.exportGraphs(prefix)
	#logging.info(greenText('Phase 4: exporting results'))
	############
	#  TODOS   #
	############
	
	# exports 
	# ???
	
