#This file generates the .txt script which can then be run by SLIM to simulate samples 


library("optparse")
# Option parsing remains the same
option_list = list(
  make_option(c("-m", "--mu"), type="numeric", default="1e-05", 
              help="mutation rate [default= 1e-06]", metavar="character"),
  make_option(c("-r", "--rho"), type="numeric", default="1e-05", 
              help="recombination rate [default= 1e-05]", metavar="character"),
  make_option(c("-N", "--Ne"), type="numeric", default="1e04", 
              help="Effective population size [default= 1e04]", metavar="character"),
  make_option(c("-M", "--samplingdepth"), type="numeric", default="800", 
              help="Sampling depth [default= 800]", metavar="character"),
  make_option(c("-i", "--rep"), type="numeric", default=1, 
              help="rep", metavar="character"),
  make_option(c("-s", "--selection"), type="character", default='neutral', 
              help="selective regime [default = 'neutral', options = 'neutral, MPL']", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Use here package for path management
library(here)

# Input args - using relative paths
filepath <- here("simulation")
mu <- opt$mu
rho <- opt$rho
Ne <- opt$Ne
M <- opt$samplingdepth
rep <- opt$rep
sel <- opt$selection

# Fixed args - using relative paths
path_to_hxb2_env <- here("HXB2.txt")  
samplingGenerations <- 5*Ne + seq(0, 600, by=50)
samplingGenerations <- round(samplingGenerations, digits=0)
trial_directory <- paste0("_rep", rep, "_rho", rho)

# Directory creation using file.path for cross-platform compatibility
if (!dir.exists(filepath)) {
  dir.create(filepath)
}

# Create directory for the particular slim run
dir.create(file.path(filepath, trial_directory))

# Rest of the script remains the same
M=300
late_string <- ""
for (i in samplingGenerations){
  late_string <- (paste0(late_string, i,  " late() { sim.outputFixedMutations(filePath=\"gen_", i, "_fixed.txt\"); } \n",
                         i, " late() { p1.outputSample(",M , ", filePath=\"gen_", i, "_seg.txt\"); } \n"))
}

# Write the slim file
slimscript <- paste0('initialize() {
 initializeSLiMOptions(nucleotideBased=T); 
 defineConstant("L", initializeAncestralNucleotides("',path_to_hxb2_env,'" ));
 initializeMutationTypeNuc("m1", 0.5, "f", 0.0);
 initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(',mu/3,'));
 initializeGenomicElement(g1, 0, L-1);
 initializeRecombinationRate(',rho,');
                        
}
1 early() {   sim.addSubpop("p1",', Ne,'); }\n \n',  
late_string
)

# Write using file.path for the output
write(slimscript, file.path(filepath, trial_directory, "slim_script.txt"))