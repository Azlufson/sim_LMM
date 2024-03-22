library(afex)
library(parallel)
library(future.apply)
library(foreach)
library(doRNG)

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b2*cond + (1|subj) + epsilon
sim_data_int <- function(n.subj = 10, n.obs = 10, b0 = 10, beta_cond = 0, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs * 2)
  cond <- rep(c(0,1), n.subj*n.obs)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs*2)
  y <- 10 + beta_cond * cond + subj_int + rnorm(length(subj), 0, 2)
  return(data.frame(subj, cond, y))
}

test_PB.fixed <- function(model, data, nsim.pb = 1000, cl = NULL) {
  return(suppressMessages(mixed(model, data = data, method = "PB", progress = FALSE, cl = NULL, args_test = list(nsim = nsim.pb, cl = cl))$anova_table$`Pr(>PB)`[1]))
}

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)



nsim.mixed <- 100 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 100


(nc <- detectCores()) # number of cores
cl <- makeCluster(rep("localhost", nc)) # make cluster

data_PB <- t(apply(grid, 1, function(x) replicate(nsim.mixed, test_PB.fixed(model, data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), nsim.pb = nsim.pb, cl = NULL))))
stopCluster(cl)
?mixed
#cl = cl, args_test = list(nsim = nsim.pb, cl = cl)
#Auslastung RAM: ~50%

#cl = NULL, args_test = list(nsim = nsim.pb, cl = cl)
#Auslastung RAM: 38% -> the same as when using no parallelization

#cl = cl, args_test = list(nsim = nsim.pb, cl = NULL)
#Auslastung RAM: 54%



###Test pbmodcomp

library(pbkrtest)
library(lme4)

m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

data <- sim_data_int()
large <- lmer(m.full, data)
small <- lmer(m.null, data)

PBmodcomp(largeModel = large, smallModel = small, nsim = 10000, cl = NULL)

(nc <- detectCores()) # number of cores
cl <- makeCluster(rep("localhost", nc)) # make cluster
PBmodcomp(largeModel = large, smallModel = small, nsim = 10000, cl = cl)
stopCluster(cl)


#cl = NULL
#Auslastung: ~37%

#find source code
methods(PBmodcomp)
pbkrtest:::PBmodcomp.merMod
getAnywhere(PBrefdist)
methods(PBrefdist)
pbkrtest:::PBrefdist.merMod
#pbkrtest:::getLRT.lmerMod
?getME.merMod #used for checking if REML == TRUE
?do_sampling

getAnywhere(do_sampling)
getAnywhere(get_refdist)
pbkrtest:::get_refdist.merMod

#trying exactLRT
data = sim_data_int()
m1 <- lmer(m.full, data = data)
m0 <- lmer(m.null, data = data)

exactLRT(m1, m0) #can only be used to test random effects

exactRLRT(m = m1, m0 = m0) #only tests random effect?


#using simulate()
#1. simulate data from fitted model
#2. fit model using new data
#3. compute LRT statistic
#4. repeat n.pb times
#5. p-value = mean(lrstat > 2.568) (cf. Faraway, p. 205)

#Funktion fürs parametric bootstrapping
#benötigt package foreach
#mixed (bz.w pbmodcomp) parallelisiert nicht
#doRNG notwendig, um richtig seed zu setzen
test_pb <- function(data, m.full, m.null, n.pb, REML = FALSE) {
  nullmod <- lmer(m.null, data = data, REML = REML)
  fullmod <- lmer(m.full, data = data, REML = REML)
  lrstat <- numeric(n.pb)
  print("1")
  lrstat <- foreach(i = 1:n.pb, .combine = "c") %dopar% {
    data$y_sim <- unlist(simulate(nullmod))
    null <- lmer(update(m.null, y_sim ~ .), data = data, REML = REML)
    alt <- lmer(update(m.full, y_sim ~ .), data = data, REML = REML)
    as.numeric(2*(logLik(alt)-logLik(null)))
  }
  return(mean(lrstat > as.numeric(2*(logLik(fullmod)-logLik(nullmod)))))
}

#test_pb(data = sim_data_int(), m.full = m.full, m.null = m.null, n.pb = 100, REML = FALSE)

n.sim <- 12
n.bs <- 100

n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

registerDoRNG(123)
start <- Sys.time()
lrt <- replicate(n.sim, test_PB(data = sim_data_int(beta_cond = 0), m.full = m.full, m.null = m.null, n.bs = n.bs, REML = FALSE))
end <- Sys.time()
end-start
parallel::stopCluster(cl = my.cluster)

beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

start <- Sys.time()
data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = FALSE))))
end <- Sys.time()
end-start
parallel::stopCluster(cl = my.cluster)


