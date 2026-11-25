seed = -1
seqfile = ../Cervinae_10inds_1out.concat.phy
treefile = ../Cervinae_10inds_1out_astral.tree
mcmcfile = mcmc_usedata3.txt
outfile = mcmc_out_usedata3.txt

ndata = 1
seqtype = 0
usedata = 3
clock = 2

model = 7
alpha = 0.5
ncatG = 5

cleandata = 0

BDparas = 0.01 0.01 0.1 M
kappa_gamma = 6 2
alpha_gamma = 1 1

rgene_gamma = 2 2000 1
sigma2_gamma = 1 10 1

finetune = 1: 0.1 0.1 0.1 0.1 0.1 0.1

print = 1
burnin = 100000
sampfreq = 1000
nsample = 100000
