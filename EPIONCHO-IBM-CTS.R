
#function rotate matrix, used in mf function 
rotate <- function(x) t(apply(x, 2, rev))

## CTS-specific function
update.dur.preg <- function(dur.preg, cur.preg, gest, death.vec, cur.age, real.max.age, DT) {
  dur.preg <- dur.preg + DT
  dur.preg[cur.preg==0 | dur.preg>gest] <- 0
  dur.preg[which(death.vec == 1)] <- 0 #set individuals which die to age 0
  dur.preg[which(cur.age >= real.max.age)] <- 0
  dur.preg[which(cur.age > 40)] <- 0
  dur.preg[which(cur.age < 16)] <- 0
  dur.preg
}

## CTS-specific function
update.preg <- function(N, p.preg, cur.age, sex, gest, dur.preg, cur.preg, death.vec, real.max.age, DT)
{
  get.preg <- rbinom(N, 1, p.preg*DT)
  get.preg[cur.age<16 | cur.age>40 | cur.preg==1 | sex==1] <- 0
  cur.preg[cur.preg==0 & get.preg==1] <- 1
  cur.preg[dur.preg>gest] <- 0
  
  cur.preg[which(death.vec == 1)] <- 0 #set individuals which die to age 0
  cur.preg[which(cur.age >= real.max.age)] <- 0
  cur.preg[which(cur.age > 40)] <- 0
  cur.preg[which(cur.age < 16)] <- 0
  
  cur.preg
}

## Contains CTS code
change.micro <- function(dat, num.comps, mf.cpt, num.mf.comps, ws, DT, time.each.comp, mu.rates.mf, fec.rates, mf.move.rate,
                         up, kap, iteration, treat.vec, give.treat, treat.start, nfix, N, do.clin.trial, start.trial)
  
{
 # N <- length(dat[,1]) 
  #indexes for fertile worms (to use in production of mf)
  fert.worms.start <-  ws + num.comps*2 
  fert.worms.end <-  (ws-1) + num.comps*3
  
  #indexes to check if there are males (males start is just 'ws')
  mal.worms.end <- (ws-1) + num.comps
  mf.mu <- rep(mu.rates.mf[mf.cpt], N)
  fert.worms <- dat[, fert.worms.start:fert.worms.end]
  
  ################################################
  ## MDA/CTS treatment function
  ################################################
  if( (give.treat == 1 & iteration >= treat.start) | (do.clin.trial==T & iteration >= start.trial) )
  {
    #which.treat <- which(is.na(treat.vec) != TRUE) #if false, the criteria for at least one treatment have been passed in the adult worm func, treat.vec[ind] is interation*DT
    ## this takes treament vec which is the most recent time since last treatment (either by CT or MDA)
    tao <- ((iteration-1)*DT) - treat.vec #tao is zero if treatment has been given at this timestep
    #print(treat.vec)
    mu.mf.prime <- ((tao + up) ^ (-kap))
    #print(mu.mf.prime)
    
    ## this sets any excexx mortality terms to 0 for individuals 
    ## not receiving treatment because of NA code
    mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0
    
    #print(sum(mu.mf.prime))
    
    mf.mu <- mf.mu + mu.mf.prime
    
  }
  ################################################
  ## end MDA treatment function
  ################################################
  
  if(mf.cpt == 1)
  {
    mp <- rep(0, N)
    
    inds.fec <- which(rowSums(dat[, ws : mal.worms.end]) > 0); mp[inds.fec] <- 1     #need to check there is more than one male
    
    k1 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, nfix + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = 0)  #fert worms and epin are vectors
    k2 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, nfix + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k1/2) 
    k3 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, nfix + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k2/2)  
    k4 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, nfix + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k3) 
    
    out <- dat[, nfix + mf.cpt] + DT/6*(k1+2*k2+2*k3+k4)
    
  }
  
  
  if(mf.cpt > 1)
  {
    k1 <- derivmf.rest(mf.in = dat[, nfix + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = 0) 
    k2 <- derivmf.rest(mf.in = dat[, nfix + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k1/2) 
    k3 <- derivmf.rest(mf.in = dat[, nfix + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k2/2) 
    k4 <- derivmf.rest(mf.in = dat[, nfix + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k3) 
    
    out <- dat[, nfix + mf.cpt] + DT/6*(k1+2*k2+2*k3+k4)
    
    #print(DT/6*(k1+2*k2+2*k3+k4))
  }  
  
  #if(out < 0) {out <- 0}
  
  return(out)
}

## No CTS code
derivmf.one <- function(fert.worms, mf.in, ep.in, mf.mort, mf.move, mp, k.in)  #fert worms and epin are vectors
{
  
  new.in <- (rotate(fert.worms)*ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each row
  new.in <- rotate(rotate(rotate(new.in)))
  new.in <- rowSums(new.in)
  
  out <- mp * new.in - mf.mort*(mf.in + k.in) - mf.move * (mf.in + k.in) 
  return(out)
}

## No CTS code
derivmf.rest <- function(mf.in, mf.mort, mf.move, mf.comp.minus.one, k.in) 
{
  out <- mf.comp.minus.one*mf.move - mf.mort*(mf.in + k.in) - mf.move * (mf.in + k.in)
  return(out)
}

#prop of L3 larvae developing into adult worms in one human, expos = total exposure for an individual
## No CTS code
delta.h <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos)
{
  out <- (delta.hz + delta.hinf * c.h * m * beta *  L3*expos) / (1 + c.h * m * beta * L3*expos)
  return(out)
}

#L1, L2, L3 dynamics
## No CTS code
calc.L1 <- function(beta, mf, mf.delay.in, expos, delta.vo, c.v, nuone, mu.v, a.v, expos.delay)
{
  delta.vv <- delta.v(delta.vo, c.v, mf, expos)
  
  #out <- (delta.vv * beta * expos *  mf)  / (mu.v + a.v * mf + nuone) 
  
  out <- (delta.vv * beta * expos *  mf)  / ((mu.v + a.v * mf) + (nuone * exp (-(4/366) * (mu.v + (a.v * mf.delay.in*expos.delay)))))
  
  #if(out < 0) {print('negative L1'); out <- 0}
  
  return(out)
}

## No CTS code
calc.L2 <- function(nuone, L1.in, mu.v, nutwo, mf, a.v)
  
{
  #out <- (nuone * L1.in) / (mu.v + nutwo) 
  
  out <- (L1.in * (nuone * exp (-(4/366) * (mu.v + (a.v * mf))))) / (mu.v + nutwo)
  
  #if(out < 0) {print('negative L2'); out <- 0}
  
  return(out)
}

## No CTS code
calc.L3 <- function(nutwo, L2.in, a.H, g, mu.v, sigma.L0)
  
{
  out <- (nutwo * L2.in) / ((a.H / g) + mu.v + sigma.L0) 
  
  #if(out < 0) {print('negative L3'); out <- 0}
  
  return(out)
}

#rate of acquisition of new infections in humans
## No CTS code∑
Wplus1.rate <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos, DT, expos.main)
  
{
  
  dh <- delta.h(delta.hz, delta.hinf, c.h, L3, m , beta, expos)
  
  out <- DT * m * beta * dh * expos.main * L3
  
  return(out)
  
}

#proportion of mf per mg developing into infective larvae within the vector
## No CTS code
delta.v <- function(delta.vo, c.v, mf, expos)
  
{
  
  out <- delta.vo / (1 + c.v * mf*expos)
  
  return(out)
  
}

#calculate number of mf in skin snip
## No CTS code
mf.per.skin.snip <- function(ss.wt, num.ss, slope.kmf, int.kMf, data, nfw.start, fw.end,  ###check vectorization 
                             mf.start, mf.end, pop.size)
  
{
  
  all.mfobs <- c()
  
  kmf <- slope.kmf * (rowSums(data[,nfw.start:fw.end])) + int.kMf #rowSums(da... sums up adult worms for all individuals giving a vector of kmfs
  
  
  mfobs <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))
  
  nans <- which(mfobs == 'NaN'); mfobs[nans] <- 0
  
  if(num.ss > 1)
    
  {
    
    tot.ss.mf <- matrix(, nrow = length(data[,1]), ncol = num.ss)
    tot.ss.mf[,1] <- mfobs
    
    for(j in 2 : (num.ss)) #could be vectorized
      
    {
      
      temp <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))
      
      nans <- which(temp == 'NaN'); temp[nans] <- 0
      
      tot.ss.mf[,j] <- temp
      
    }
    
    mfobs <- rowSums(tot.ss.mf)  
    
  } 
  
  mfobs <- mfobs / (ss.wt * num.ss) 
  
  list(mean(mfobs), mfobs)
  
}

#calculate CMFL from skin snip output
## No CTS code
calc.CMFL <- function(ss.mf.pi, data, ss.wt)
  
{
  ind.ov.20 <-  which(data[,2] > 20)
  mf.o.t <- log((ss.mf.pi[ind.ov.20] * ss.wt) + 1)
  CMFL <- exp(sum(mf.o.t)/length(ind.ov.20)) - 1  
  
  return(CMFL)
}

#change in the number of adult worms in an indiviudal host, starts with males then moves to females (new female worms are calculated in male 'section')
## Contains CTS components
change.worm.per.ind<- function(delta.hz, delta.hinf, c.h, L3, m , beta, compartment, total.dat, num.comps,
                               w.f.l.c, lambda.zero, omeg, expos, ws, DT, mort.rates, mort.rates.macro, time.each.comp, new.worms.m, new.worms.nf.fo,
                               rn, lam.m, phi, treat.stop, iteration, treat.int, treat.prob, cum.infer, treat.vec, give.treat, 
                               treat.start, N, ageindi, compindi, cohortindi, do.clin.trial, MOM, macro.t, start.trial, stop.trial, trial.int, stop.treat.trial, MDA.vec,  CT.vec)
  
{
  
  ############################################################################################################
  ### CTS code chunck -  swapping background macrofilaricidal rates for increased rates in treated individuals
  ############################################################################################################
  
  
  #N <- length(treat.vec)
  
  #print(N)
  
  lambda.zero.in <- rep(lambda.zero * DT, N) #loss of fertility
  omeg <- rep(omeg * DT, N) #becoming fertile
  
  #male worms
  
  cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment
  
  cur.Wm <- total.dat[, cl] #take current number of worms from matrix
  
  #female worms
  
  clnf <- (ws - 1) + num.comps + compartment #column for infertile females, num.comps skips over males
  
  clf <- (ws - 1) + 2*num.comps + compartment #column for fertile females, 2*num.comps skips over males and infertile females
  
  cur.Wm.nf <- total.dat[, clnf] #take current worm number from matrix infertile females 
  
  cur.Wm.f <- total.dat[, clf] #take current worm number from matrix fertile females
  
  mort.fems <- mort.mals <- rep(mort.rates[compartment], N)
  mort.macro <- rep(mort.rates.macro[compartment], N)
  
  #####################################################  
  ######### 
  #treatment
  #########
  
  #check if a complier, older than five, time for treatment, still within treatment timeframe
  #approach assumes individuals which are moved from fertile to non fertile class due to treatment re enter fertile class at standard rate (bit wierd)
  
  ########################################################
  ### MDA/CTS code chunck -  could be written as a function
  #######################################################
  if(give.treat == 1 & iteration >= treat.start) ####CHECK 
  {
    # identify treatment times
    if((((iteration-1)*DT) %% treat.int < DT) & iteration <= treat.stop) #is it a treatment time T/F
    {
      
      inds.to.treat <- which(total.dat[,compindi]  == 0 & 
                               total.dat[,ageindi] > 5 & rn <= treat.prob) #find individuals which are compliers, older than 5 and under the coverage
      MDA.vec[inds.to.treat]  <-  (iteration-1) * DT #alter vector storing time of treatment 
      if(iteration > treat.start) {mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer)}#; print('cum.infer')}
      #mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer) #alter mortality for cumulative effects of IVM
      #mort.mals[inds.to.treat] <- mort.mals[inds.to.treat] + (cum.infer)
    }
  }
  if(do.clin.trial==TRUE & iteration >= start.trial) {
      # is it a trial treatment time?
      # this code means that treatment starts 1 step size after iteration=start.trial
      if(((iteration -1 - start.trial) %% trial.int < DT) & iteration <= stop.treat.trial) 
      {
      print(iteration)
        # cohort 1 receives ivermectin 
        inds.to.treat <- which(total.dat[,cohortindi]  == 1) #find individuals assigned to control group cohort 1
        if(iteration > start.trial) {mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer)}#; #print('cum.infer')}
       # mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer) #alter mortality for cumulative effects of IVM
       # mort.mals[inds.to.treat] <- mort.mals[inds.to.treat] + (cum.infer)
        
        ## record time since lastreatment for cohort 1 and 2 - 2 required for duration of macro activty
        inds.to.treat <- which(total.dat[,cohortindi]  == 1 | total.dat[,cohortindi]  == 2) #find individuals assigned to control group cohort 1
        ## or tested group 2
        CT.vec[inds.to.treat]  <-  (iteration-1) * DT #alter time since treatment for concomitant microfilaricidal activity 
      
      }
    
    tmp.tsince <- (iteration-1)*DT - CT.vec
  #  print(tmp.tsince)
    inds.to.macro <- which(total.dat[,cohortindi]  == 2 & tmp.tsince <= macro.t ) #find individuals for which macro effect should be applied
  #  print(inds.to.macro)
    mort.fems[inds.to.macro] <- mort.macro[inds.to.macro] #macrofilaricidal activity without cumulative adjustment
    mort.mals[inds.to.macro] <- mort.macro[inds.to.macro] 
    
   # print(paste("treated mortality rate", mean(mort.fems[inds.to.macro])))
  #  print(paste("control mortality rate", mean(mort.fems[-c(inds.to.macro)])))
  }
  
  ## set a temoraary CT.vec for MOM activity so that microfilaricidal effects are not modelled
  CT.vec.tmp <- CT.vec
  if (MOM) {
    inds.to.treat <- which(total.dat[,cohortindi]  == 2) 
    CT.vec.tmp[inds.to.treat] <- NA
  }
  
  # treat.vec for output to increase mf mortality rates
  # take the most recent time
  treat.vec <- pmax(MDA.vec, CT.vec.tmp, na.rm=T)
  # treat vec to adjust female fertility rates being unaffected by macrofilaricide
  CT.vec.fert <- CT.vec
  CT.vec.fert[which(total.dat[,cohortindi]  == 2)] <- NA
  treat.vec.fert <- pmax(MDA.vec, CT.vec.fert, na.rm=T)
  
  ############################################################
  ## end MDA/CTS code chunck
  ############################################################
  
  ## simulate death and ageing of males
  worm.dead.males <- rbinom(N, cur.Wm, mort.mals)
  worm.loss.males <- rbinom(N, (cur.Wm - worm.dead.males), rep((DT / time.each.comp), N))
  
  if(compartment == 1)
    
  {
    
    male.tot.worms <- cur.Wm + new.worms.m - worm.loss.males - worm.dead.males
  }
  
  if(compartment > 1)
    
  {
    male.tot.worms <- cur.Wm + w.f.l.c[[2]] - worm.loss.males - worm.dead.males
  }
    
  ##simulate fertility movement of females
    tao <- ((iteration-1)*DT) - treat.vec.fert #time since last treatment affecting fertility
    
   # print(mean(tao))
  
    lam.m.temp <- rep(0, N); lam.m.temp[which(is.na(treat.vec.fert) != TRUE)] <- lam.m
    
    f.to.nf.rate <- DT * (lam.m.temp * exp(-phi * tao)) 
    
    f.to.nf.rate[which(is.na(treat.vec.fert) == TRUE)] <- 0 #these entries in f.to.nf.rate will be NA, lambda.zero.in cannot be NA
    
    lambda.zero.in <- lambda.zero.in + f.to.nf.rate #update 'standard' fertile to non fertile rate to account for treatment 
    

  
  #.fi = 'from inside': worms moving from a fertile or infertile compartment
  #.fo = 'from outside': completely new adult worms 
  
  ##simulate death and aging of females
  worm.dead.nf <- rbinom(N, cur.Wm.nf, mort.fems) #movement to next compartment
  
  worm.dead.f <- rbinom(N, cur.Wm.f, mort.fems)
  
  worm.loss.age.nf <- rbinom(N, (cur.Wm.nf - worm.dead.nf), rep((DT / time.each.comp), N))
  
  worm.loss.age.f <- rbinom(N, (cur.Wm.f - worm.dead.f), rep((DT / time.each.comp), N))
  
  
  #calculate worms moving between fertile and non fertile and deaths and aging 
  
  #females from fertile to infertile
  
  new.worms.nf.fi <- rep(0, N)
  
  trans.fc <- which((cur.Wm.f - worm.dead.f - worm.loss.age.f) > 0)
  #print(trans.fc)
  
  if(length(trans.fc) > 0)
  {
    new.worms.nf.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.f[trans.fc] - worm.dead.f[trans.fc] - worm.loss.age.f[trans.fc]), lambda.zero.in[trans.fc]) 
  }
  
  #else {new.worms.nf.fi <- 0}
  
  #females worms from infertile to fertile #this happens independent of males, but production of mf depends on males
  
  new.worms.f.fi <- rep(0, N)
  
  trans.fc <-  which((cur.Wm.nf - worm.dead.nf - worm.loss.age.nf) > 0)
  if(length(trans.fc) > 0)
  {
    new.worms.f.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.nf[trans.fc] - worm.dead.nf[trans.fc] - worm.loss.age.nf[trans.fc]), omeg[trans.fc])#females moving from infertile to fertile
  }
  
  if(compartment == 1)
    
  {
    nf.out <- cur.Wm.nf + new.worms.nf.fo + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi - worm.dead.nf#final number of infertile worms
    
    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi - worm.dead.f#final number of fertile worms
  }       
  
  if(compartment > 1)
    
  {
    nf.out <- cur.Wm.nf + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi + w.f.l.c[[5]] - worm.dead.nf#w.f.l.c = worms from previous compartment
    
    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi + w.f.l.c[[6]] - worm.dead.f
  }   
  
  # if(male.tot.worms < 0) {male.tot.worms <- 0}
  # if(nf.out < 0) {nf.out <- 0}
  # if(f.out < 0) {f.out <- 0}
  
  list(male.tot.worms,
       worm.loss.males,
       nf.out,
       f.out,
       worm.loss.age.nf,
       worm.loss.age.f, treat.vec, MDA.vec, CT.vec)  
}

weibull.mortality <- function(DT, par1, par2, age.cats)
  
{
  ## reparameterized - default for adult worms is 
  ## take midpoint of age categories for hazard rates
  d <- diff(age.cats)[1]/2
  out <- DT  * (par1 ^ (par2)) * par2 * ( (age.cats+d) ^ (par2-1) )
  
  return(pmin(out, 1))
  
  #return(out)
}

prevalence.for.age <- function(age, ss.in, main.dat)
  
{
  
  inds <- which(main.dat[,2] >= age)
  
  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds) 
  
  return(out)
}

run.mod <- function(ABR=1000,delta.hz  = 0.1864987,delta.hinf = 0.002772749 , c.h = 0.004900419,
                        m.exp = 1.08, f.exp = 0.9, age.exp.m = 0.007, age.exp.f = -0.023, mu.v = 26, int.mf = 0, sigma.L0 = 52,
                        a.H = 0.8, g = 0.0096, h=0.63, a.v = 0.39, num.comps.worm = 20, real.max.age = 80, N = 100, mean.age = 50,
                        int.L3 = 0.03, int.L2 = 0.03, int.L1 = 0.03, lambda.zero = 0.33, omeg = 0.59,
                        delta.vo = 0.0166, c.v = 0.0205, num.mf.comps = 20, DT = 2/366, int.worms=1, ss.wt = 2, num.ss = 2,
                        slope.kmf = 0.0478, int.kMf = 0.313, sex.rat = 0.5, nuone = 201.6189, nutwo = 207.7384, 
                        time.each.comp.worms = 1, time.each.comp.mf = 0.125, mu.w1 = 0.0975, mu.w2 = 4.25, fec.w.1 = 70, fec.w.2 = 0.72,
                        mu.mf1 = 0.95, mu.mf2 = 1.6,  sero.sensitivity.in = 0.8, mf.move.rate = 8.133333,
                        type.sero = 0, l3.delay = 10, lam.m = 32.4, phi = 19.6, stop.time=50, treat.prob = 0.75,
                        cum.infer= 0.345, up = 0.0096, kap = 1.25, give.treat = 1, treat.start = 15, n.trt=5, treat.int=1, 
                        every = 2/366, pnc = 0.05, gest=2, p.preg=0.05, ke=0.3, age.prev.in=5,
                        ex.dd.int=0, get.equib=0, input.var = 0, input.eq=NA, do.clin.trial=1,
                        clin.trial=list(MOM=1, macro.t=3/12, macro.eff=0.9, start.trial=1, stop.trial=2, n.trt.trial=NA, trial.int=NA, min.age=NA, 
                                        target.ss=NA, min.mf=5, drop.out=0.1), print.progress=1)
{ 
  ## derived variables
  beta = h/g
  m = ABR*(g / h)
  
  dt.days <- DT*366
  time.its = round(stop.time / DT)
  
  it.thresh=time.its+1 ## this needs to be recoded to allow serology to be measured ##
  
  treat.start = round(treat.start / (DT)); treat.stop = treat.start + round(n.trt*(treat.int/(DT )))
  
  #gam.dis <- (ABR ^ ex) * cor + int  #theta / k for gamma
  gam.dis <- ke
  
  ###################################################################################
  ## CTS-specific code - extract list element parameters & define indicator variables
  ###################################################################################
  ## do not do clinical trial if equilibrium inputs are included (this could be modified)
  ## to be conditional on specific equilibrium inputs
  if (do.clin.trial) {
    MOM <- clin.trial$MOM
    macro.t <- clin.trial$macro.t
    macro.eff <- clin.trial$macro.eff
    start.trial = round(clin.trial$start.trial / DT)
    trial.int <- clin.trial$trial.int/DT
    stop.treat.trial = start.trial + round(clin.trial$n.trt.trial*(trial.int))
    stop.trial = round(clin.trial$stop.trial / DT)
    target.ss <- clin.trial$target.ss
    min.mf <- clin.trial$min.mf
    min.age <- clin.trial$min.age
    ## per year drop out rate
    drop.out.rate <- -log(1-clin.trial$drop.out)
    ## probability of drop out per time step
    prob.drop <- drop.out.rate*DT
    compindi <- 1
    ageindi <- 2
    sexindi <- 3
    pregindi <- 4
    durpregindi <- 5
    cohortindi <- 6
    L1indi <- 7
    L2indi <- 8
    L3indi <- 9
    nfix <- L3indi ## nfix indicates single column variables
    N.cohorts <- 2
  } else {
    MOM <- clin.trial$MOM
    macro.t <- clin.trial$macro.t
    macro.eff <- clin.trial$macro.eff
    start.trial = round(stop.time/DT) + 1
    trial.int <- clin.trial$trial.int/DT
    stop.treat.trial = start.trial + round(clin.trial$n.trt.trial*(trial.int))
    stop.trial = stop.treat.trial + 1
    target.ss <- clin.trial$target.ss
    min.mf <- clin.trial$min.mf
    min.age <- clin.trial$min.age
    ## per year drop out rate
    drop.out.rate <- -log(1-clin.trial$drop.out)
    ## probability of drop out per time step
    prob.drop <- drop.out.rate*DT
    compindi <- 1
    ageindi <- 2
    sexindi <- 3
    pregindi <- 4
    durpregindi <- 5
    cohortindi <- 6
    L1indi <- 7
    L2indi <- 8
    L3indi <- 9
    nfix <- L3indi ## nfix indicates single column variables
    N.cohorts <- 2
  }
  ############################################################
  ## end of CTS code chunk
  ###########################################################
  
  ## if running to equilibrium, turn off treatment and trial starts
  if (get.equib==1) {
    give.treat <- 0
    treat.start <- round(stop.time / DT) + 10
    start.trial <- round(stop.time / DT) + 10
    stop.trial <- start.trial + 10
  }
  
  
  #columns to set to zero when an individual dies
  
  cols.to.zero <- seq(from = 1, to = (nfix + num.mf.comps + 3*num.comps.worm))  #should this include 1 for treatment?
  cols.to.zero <- cols.to.zero[-c(1,L2indi, L3indi)] #compliance, L2 and L3 do not become zero when an individual dies
  
  #used to perform operations on different worm and mf compartments 
  tot.worms <- num.comps.worm*3
  num.cols <- nfix + num.mf.comps + tot.worms 
  worms.start <- (nfix+1) + num.mf.comps
  
  
  nfw.start <- (nfix+1) + num.mf.comps + num.comps.worm #start of infertile worms
  fw.end <- num.cols #end of fertile worms 
  mf.start <- (nfix+1)
  mf.end <- nfix + num.mf.comps
  
  #list for sero prevalence by age group
  if(type.sero == 1) {sprev.by.group.m <- list(length=12); mfprev.by.group.m <- list(length=12); sprev.by.group.f <- list(length=12); mfprev.by.group.f <- list(length=12)} 
  else {sprev.by.group.m <- 0; mfprev.by.group.m <- 0; sprev.by.group.f <- 0; mfprev.by.group.f <- 0}
  
  
  #age dependent mortality and fecundity rates
  
  age.cats <- seq(0, 20, length = num.comps.worm)
  
  ## mortality rates for each compartment
  mort.rates.worms <- weibull.mortality(DT = DT, par1 = mu.w1, par2 = mu.w2, age.cats = age.cats)
  mort.rates.worms.macro <- mort.rates.worms 
  
  if (do.clin.trial) {
   # mu.w1.macro <- exp( 1/mu.w2*log( -log(1-macro.eff)) - log(macro.t)  )
    mu.w1.macro <- min(DT*(-log(1-macro.eff)/(macro.t)), 1)
    mort.rates.worms.macro <- weibull.mortality(DT = DT, par1 = mu.w1, par2 = mu.w2, age.cats = age.cats) + mu.w1.macro
    #print(mort.rates.worms.macro)
  }
  
  #mort.rates.worms <- rep(0.1*DT, num.comps.worm) #to test no senescence 
  
  fec.rates.worms <- 1.158305 * fec.w.1 / (fec.w.1 + (fec.w.2 ^ -age.cats) - 1) #no DT - Rk4
  
  #fec.rates.worms <- rep(1.15, num.comps.worm) #to test no senescence 
  
  age.cats.mf <- seq(0, 2.5, length = num.mf.comps)
  
  mort.rates.mf <- weibull.mortality(DT = 1, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf)
  
  #mort.rates.mf <- rep(1.2, num.mf.comps) #to test no senescence 
  
  ###########################################################
  ## CTS-specific code - initialise pegnancy & output vectors 
  ###########################################################
  if (do.clin.trial) {
    dur.preg <- c()
    cur.preg <- c()
    dur.preg <- rep(0, N)
    cur.preg <- rep(0, N)
    prev.by.ss.cohort <- mn.mf.by.ss.cohort <- sd.prev.by.ss.cohort <-
      sd.mf.by.ss.cohort <- n.per.cohort <- matrix(0, nrow=time.its+1, ncol=N.cohorts)
    prev.diff <- mf.diff <- rep(0, length=time.its+1)
  }
  ###########################################################
  ## end of CTS code chunk 
  ###########################################################
  
  if(input.var == 0) #are we inputting no equilibirum condition?
    
  {
    
    sero.prev.vec <- 0
    #sero sensitivity 
    
    num <- ceiling(sero.sensitivity.in * N)
    ind.sensit <- sample(seq(1, N), num)
    
    sensitivity <- rep(0, N)
    sensitivity[ind.sensit] <-  1 #infected individuals either test positive or not for life
    
  
    ## assign sex
    sex <- rbinom(N, 1, sex.rat)
    
    
    #create inital age distribution & pregnancy status and simulate stable age distribution
    
    cur.age <- rep(0, N)

    
    for(i in 1 : 75000) #if at equilibrium you saved the age at which inds die and simulated further, you should get an exponential distribution
    {
      cur.age <- cur.age + DT
      
      death.vec <- rbinom(N, 1, (1/mean.age) * DT) 
      
      cur.age[which(death.vec == 1)] <- 0 #set individuals which die to age 0
      cur.age[which(cur.age >= real.max.age)] <- 0
      
      #########################################################
      ## CTS-specific code - run pregnancy model to equilibrium
      #########################################################
      if(do.clin.trial) {
        dur.preg <- update.dur.preg(dur.preg=dur.preg, cur.preg = cur.preg, gest=gest,
                                    death.vec=death.vec,cur.age=cur.age, real.max.age=real.max.age, DT=DT)
        cur.preg <- update.preg(N=N, p.preg=p.preg, cur.age=cur.age, sex=sex, gest=gest, 
                                dur.preg=dur.preg, cur.preg=cur.preg, death.vec=death.vec, real.max.age=real.max.age, DT=DT)
      }
      ###########################################################
      ## end of CTS code chunk
      ###########################################################
      
    }
    
    
    #create list to store matrices (where required) and mean number of parasites per person
    
    all.mats <- vector("list", length = time.its + 1) 
    
    #create initial exposure vector
    
    ex.vec <-rgamma(N, gam.dis, gam.dis)
    
    
    ###############################################
    #matrix for delay in l3 establishment in humans 
    num.delay.cols <- l3.delay * (28 / dt.days) 
    l.extras <- matrix(0, ncol= num.delay.cols, nrow= N)
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time
    ################################################
    
    
    
    ###############################################
    #matrix for l1 delay
    
    # num.lmd.cols <- 4 / dt.days
    # l1.delay <- matrix(int.L1, ncol= num.lmd.cols, nrow= N)
    # inds.d.mats <- seq(2,(length(l1.delay[1,])))
    
    
    l1.delay <- rep(int.L3, N)
    
    ###############################################
    
    ###############################################
    #matrix for tracking mf for l1 delay
    
    
    num.mfd.cols <- 2
    # 
    mf.delay <- matrix(int.mf, ncol= num.mfd.cols, nrow= N)
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))
    
    #mf.delay <- rep(int.mf, N)
    ###############################################
    #matrix for exposure for delay
    num.exp.cols <- 2
    exposure.delay <- matrix(ex.vec, ncol= num.exp.cols, nrow= N)
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))  
    
    #matrix for first time step
    all.mats[[1]] <- matrix(, nrow=N, ncol=num.cols)
    
    all.mats[[1]][,  (worms.start) : num.cols] <- int.worms
    
    all.mats[[1]][, L1indi] <- int.L1
    
    all.mats[[1]][, L2indi] <- int.L2
    
    all.mats[[1]][, L3indi] <- int.L3
    
    all.mats[[1]][, (nfix+1) : (nfix + 1 + (num.mf.comps-1))] <- int.mf
    
    all.mats[[1]][,compindi] <- rep(0, N) #column used during treatment
    all.mats[[1]][,ageindi] <- cur.age
    
    #############################################################################
    ## CTS-specific code - store pregnancy status
    #############################################################################
    if (do.clin.trial) {
      all.mats[[1]][,pregindi] <- cur.preg
      all.mats[[1]][,durpregindi] <- dur.preg
      all.mats[[1]][,cohortindi] <- rep(0, N)
    }
    ##############################################################################
    ## end of CTS code chunk
    ##############################################################################
    
    #store sex to humans 
    
    all.mats[[1]][,sexindi] <- sex
    

    #mf per snip and CMFL
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats[[1]], nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    mean.mf.per.snip <- temp[[1]]
    CMFL<- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats[[1]], ss.wt)
    
    #inds.five <- which(all.mats[[1]][,ageindi] >= 5)
    
    #prev.by.ss <- length(which(temp[[2]][inds.five] > 0)) / length(inds.five) 
    
    prev.by.ss <- prevalence.for.age(age = age.prev.in, ss.in = temp, main.dat = all.mats[[1]])
    
    
    sero.prev <- c()
    sero.prev.10 <- c()
    
    test.prev <- c()
    
    L3 <- mean(all.mats[[1]][, L3indi])
    
    L2 <- mean(all.mats[[1]][, L2indi])
    
    L1 <- mean(all.mats[[1]][, L1indi])
    
  }
  
  ###############################################################################
  ## end of part which is for equilibirum only, or treatment requiring equilibrium
  ################################################################################

  else{
    
    all.mats <- list(input.eq[[1]])
    ex.vec <- input.eq[[2]]
    sensitivity <- input.eq[[3]]
    sero.prev.vec <- input.eq[[4]]
    exposure.delay <- input.eq[[9]] 
    
    #num.exp.cols <- 2 #THIS NEEDS TO MADE FLEXIBLE WHEN DELAY SITUATION IS 'SOLVED'
    inds.exp.mats <- seq(2,(length(exposure.delay[1,]))) 
    
    l.extras <- input.eq[[5]]
    inds.l.mat <- seq(2,(length(l.extras[1,])))
    
    l1.delay <- input.eq[[6]]
    #num.mfd.cols <- 2 #THIS NEEDS TO MADE FLEXIBLE WHEN DELAY SITUATION IS 'SOLVED'
    mf.delay <- input.eq[[7]]
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))
    
    
    L3 <- mean(all.mats[[1]][, L3indi]) # put indicies in
    
    L2 <- mean(all.mats[[1]][, L2indi])
    
    L1 <- mean(all.mats[[1]][, L1indi])
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats[[1]], nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    mean.mf.per.snip <- temp[[1]]
    CMFL<- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats[[1]], ss.wt)
    
    inds.five <- which(all.mats[[1]][,ageindi] >= 5) ## check indicies
    
    prev.by.ss <- length(which(temp[[2]][inds.five] > 0)) / length(inds.five) 
    
    sero.prev <- c()
    sero.prev.10 <- c()

    if(type.sero == 1)
    {
      sero.prev.vec <- input.eq[[4]]
      sero.prev <- sum(sero.prev.vec) / N
      sero.prev.10 <- sum(sero.prev.vec[which(all.mats[[1]][,ageindi] < 10)]) / length(which(all.mats[[1]][,ageindi] < 10))
    }
    
    treat.vec <- input.eq[[8]]#for time since treatment calculations 
    
    #if(input.var == 'treatment') {iter.last.treat <- input.eq[[10]]}
    
  }
 
  nw.rate <- 1
  
  non.comp <- ceiling(N * pnc)
  out.comp <- rep(0, N)
  s.comp <- sample(N, non.comp)
  out.comp[s.comp] <- 1
  all.mats[[1]][,compindi] <- out.comp
  
  treat.vec <- MDA.vec <- CT.vec <- rep(NA, N) #for time since treatment calculations 
  
  #for debugging
  #l2par <- c()
  #l2par2 <- c()
  
  #################################################################
  ### time loop
  #################################################################
  for(i in 1:time.its) 
    
  {
    ## print information on where the simulation is
    if (print.progress==1) {
    r.time <- round(i * DT, digits = 2)
    if (r.time%%1==0) {
      print(paste(r.time, 'yrs', sep = ' '))
    }
    }
    
    #stores means from previous steps in an attempt to reduce 
    #the amount of memory used. This is a bit sloppy and needs to be replaced with a ring buffer (see ritch fitzjohn 'ring' package)
    
    if(i > 1) 
    {
      means <- vector(length=2)
      
      L3 <- mean(all.mats[[i-1]][, L3indi])
      
      if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i-1]][, (nfix+1)])))} 
      else means[1] <- mean(as.numeric(rowSums(all.mats[[i-1]][, (nfix+1):(worms.start-1)]))) #means for mf
      
      means[2] <- mean(as.numeric(rowSums(all.mats[[i-1]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )]))) #means for adult worms
      
      #############################################################################
      ## CTS-specific code - caclulate mean number of adult worms
      #############################################################################
      if (do.clin.trial) {
        means.cohort <- vector(length=2)
        means.cohort[1] <- means.cohort[2] <- 0
        if (i>=(start.trial) & i<=(stop.trial+1)) { ## removed +1
        indi.tmp <- which(all.mats[[i-1]][,cohortindi]==1)
        if (length(indi.tmp)>0) {
        means.cohort[1] <- mean(as.numeric(rowSums(all.mats[[i-1]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ),drop=F])))
        } else {means.cohort[1]<-0}
        indi.tmp <- which(all.mats[[i-1]][,cohortindi]==2)
        if (length(indi.tmp)>0) {
        means.cohort[2] <- mean(as.numeric(rowSums(all.mats[[i-1]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ),drop=F])))
        } else {means.cohort[2]<-0}
        }
        all.mats[[i-1]] <- c(means, L3, means.cohort) #replace matrix at previous timestep with means for worms and L3
        } else {
          all.mats[[i-1]] <- c(means, L3) #replace matrix at previous timestep with means for worms and L3
        }
      ##############################################################################
      ## end of CTS-specific code chunck
      ##############################################################################
      
    }
    
    all.mats.cur <- all.mats[[i]] #create temporary matrix for t rather than using list element to increase speed
    all.mats.temp <- all.mats.cur # all.mats.temp is the matrix to be updated 
    
    cur.age <- all.mats.cur[,ageindi]
    sex <- all.mats.cur[,sexindi]
    
    #sex and age dependent exposure
    
    mls <- which(all.mats.cur[,sexindi] == 1)
    fmls <- which(all.mats.cur[,sexindi] == 0)
    
    s.a.exp <- c()
    
    s.a.exp[mls] <- m.exp * exp(-age.exp.m * (all.mats.cur[mls, ageindi]))
    gam.m <- 1 / mean(s.a.exp[mls]) #normalize so mean = 1
    s.a.exp[mls] <- s.a.exp[mls] * gam.m
    
    s.a.exp[fmls] <- f.exp * exp(-age.exp.f * (all.mats.cur[fmls, ageindi]))
    gam.f <- 1 / mean(s.a.exp[fmls]) #normalize so mean = 1
    s.a.exp[fmls] <- s.a.exp[fmls] * gam.f
    
    ex.vec <- ex.vec * (1 / mean(ex.vec)) #normalize so mean = 1
    
    tot.ex.ai <- s.a.exp * ex.vec
    tot.ex.ai <- tot.ex.ai * (1 / mean(tot.ex.ai)) #normalize so mean = 1
    
    #increase age (for next time step)
    
    all.mats.temp[,ageindi] <- (all.mats.cur[,ageindi]) + DT #increase age for all individuals
    
    death.vec <- rbinom(N, 1, (1/mean.age) * DT) #select individuals to die
    
    to.die <- which(death.vec == 1)
    
    at.ab.max <- which(all.mats.temp[,ageindi] >= real.max.age)
    
    to.die <- c(to.die, at.ab.max)
    
    to.die <- unique(to.die) #may have repeated indivudals i.e selected by binom and >80
    
    
    ##############################################################
    ## CTS-specific code - update pregnancy status & select cohorts
    ###############################################################
    if (do.clin.trial) {
      dur.preg <- (all.mats.cur[,durpregindi])
      cur.preg <- all.mats.cur[,pregindi]
      dur.preg <- update.dur.preg(dur.preg=dur.preg, cur.preg=cur.preg, gest=gest,
                                death.vec=death.vec, cur.age=cur.age, real.max.age=real.max.age,
                                DT=DT)
      all.mats.temp[,pregindi] <- update.preg(N=N, p.preg=p.preg, cur.age=cur.age, sex=sex, gest=gest, 
                                            dur.preg=dur.preg, cur.preg=cur.preg, death.vec=death.vec, 
                                            real.max.age=real.max.age, DT=DT) ## update preganancy status of all individuals
      ## select the cohort two steps before the start of the trial
      if (i==(start.trial-2)) { # 2 steps before the trials 
      eligible <- rep(0,N)
      ## eligibility based on age, pregnancy status, being positive for mf
      ## and not being a non-complier!
      eligible[which(all.mats.cur[,ageindi]>min.age & 
                       all.mats.cur[,pregindi]!=1 &
                       temp[[2]]>min.mf & all.mats.cur[,compindi]==0)] <- 1
      N.eligible <- sum(eligible)
      actual.ss <- floor(N.eligible/N.cohorts) 
      applied.ss <- min(c(target.ss, actual.ss))
      ## define cohorts by dividing eligibles up untill target is met
      ## 
      indicies <- which(eligible>0)
   
      if (actual.ss>target.ss) {
        ## two cohorts
        if (length(indicies)>0) {
        indicies <- sample( indicies,target.ss*N.cohorts)
        } else {indicies<-0}
      }
      c1 <- sample(indicies, applied.ss)
      c2 <- indicies[indicies%in%c1==F]
      all.mats.temp[c1,cohortindi] <- 1 ## changed to temp
      all.mats.temp[c2,cohortindi] <- 2 ## changed to temp
      
      n.per.cohort[i+1, ] <- c(sum(all.mats.temp[,cohortindi]==1), 
                               sum(all.mats.temp[,cohortindi]==2))
      ###############test code#######
      prev.by.ss.cohort[i+1,] <- c(length(which(temp[[2]][c1] > 0)) / length(c1),
                                   length(which(temp[[2]][c2] > 0)) / length(c2))
      prev.diff[i+1] <- prev.by.ss.cohort[i, 1] - prev.by.ss.cohort[i, 2]
      ## infected or not output
      tmp <- as.numeric(temp[[2]]>0)
      # prev.by.ss.cohort[i+1,] <- c( mean(tmp[c1]), mean(tmp[c2]) )  
      sd.prev.by.ss.cohort[i+1,] <- c( sd(tmp[c1]), sd(tmp[c2]) )  
      mn.mf.by.ss.cohort[i+1,] <- c( mean(temp[[2]][c1]), mean(temp[[2]][c2]) )
      mf.diff[i+1] <- mn.mf.by.ss.cohort[i, 1] - mn.mf.by.ss.cohort[i, 2]
      sd.mf.by.ss.cohort[i+1,] <- c( sd(temp[[2]][c1]), sd(temp[[2]][c2]) )
      ############################
      
      
      
      } 
    }
    ##############################################################
    ## end of CTS code chunk
    ###############################################################
    
    ##################
    #delay calculations 
    ##################
    
    new.worms.m <- c()
    new.worms.nf <- c()
    
    new.worms.m <- rbinom(N, size = l.extras[,length(l.extras[1,])], prob = 0.5) #draw males and females from last column of delay matrix
    new.worms.nf <- l.extras[,length(l.extras[1,])] - new.worms.m
    
    #move individuals along
    
    l.extras[,inds.l.mat] <- l.extras[,(inds.l.mat-1)]
    
    #new establishing L3 vectorized for all individuals
    
    if (ex.dd.int == 1) {in.e <- tot.ex.ai} else {in.e <- rep(1, N)}
    
    nw.rate.tminus1 <- nw.rate
    
    L3.in <- mean(all.mats.cur[, L3indi])
    
    nw.rate <- Wplus1.rate(delta.hz, delta.hinf, c.h, L3.in, m ,
                           beta, expos = in.e, DT, expos.main = tot.ex.ai)
    
    new.worms <- rpois(N, nw.rate) #total new establishing L3 for each individual 
    
    l.extras[,1] <- new.worms
    
    rann <- runif(N, 0, 1) #random number for each individual to see if treatment is given (depending on compliance)
    
    if(i == 1) {treat.vec.in <- treat.vec 
    MDA.vec.in = MDA.vec
    CT.vec.in = CT.vec} #debugging ??
    
    for(k in 1 : num.comps.worm) #go through each adult worm compartment
      
    {
      
      if(k == 1) {from.last <- rep(0, N)} #create vector for worms coming from previous compartment (needs to be 0 when k ==1)
      
      res <- change.worm.per.ind(delta.hz=delta.hz, delta.hinf=delta.hinf, c.h=c.h, L3 = L3.in, m=m , beta=beta, compartment = k, 
                                 total.dat = all.mats.cur, num.comps = num.comps.worm,
                                 w.f.l.c = from.last, lambda.zero=lambda.zero, omeg=omeg, expos = tot.ex.ai, 
                                 ws = worms.start, DT = DT, mort.rates = mort.rates.worms, mort.rates.macro=mort.rates.worms.macro, time.each.comp = time.each.comp.worms, new.worms.m = new.worms.m, 
                                 new.worms.nf.fo = new.worms.nf, rn = rann, lam.m=lam.m, phi=phi, treat.stop=treat.stop, iteration = i, 
                                 treat.int=treat.int, treat.prob=treat.prob, cum.infer=cum.infer, treat.vec = treat.vec.in, 
                                 give.treat=give.treat, treat.start=treat.start,N=N, ageindi=ageindi, compindi=compindi, 
                                 cohortindi=cohortindi, do.clin.trial=do.clin.trial, MOM=MOM, macro.t=macro.t ,start.trial=start.trial,
                                 stop.trial=stop.trial, trial.int=trial.int,  stop.treat.trial=stop.treat.trial, MDA.vec=MDA.vec.in, CT.vec=CT.vec.in)
      
      from.last <- res #assign output to use at next iteration, indexes 2, 5, 6 (worms moving through compartments)
      
      #update male worms in matrix for compartment k
      
      all.mats.temp[, (nfix + num.mf.comps + k)] <- res[[1]]
      
      #update females worms in matrix
      
      all.mats.temp[, (nfix + num.mf.comps + num.comps.worm + k)] <- res[[3]] #infertile, num.comps.worm skips over males
      all.mats.temp[, (nfix + num.mf.comps + 2*num.comps.worm + k)] <- res[[4]] #fertile, num.comps.worm skips over males and infertile females
      
      #print(head(all.mats.temp))
      
    }
    
    ## update treatment vector depdening on last MDA or CT treatment
    if( (give.treat == 1 & i >= treat.start) | (do.clin.trial==T & i >= start.trial) ) {
      treat.vec.in <- res[[7]]
      MDA.vec.in <- res[[8]]
      CT.vec.in <- res[[9]]
     #print(mean(treat.vec.in, na.rm=T))
     } # update treatment vector indicating time 
    ## since last treatment (by clinical trial or MDA)
    
    for(mf.c in 1 : num.mf.comps)   
      
    {
      
      res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c, 
                             num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf, 
                             mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate=mf.move.rate, up=up, kap=kap, iteration = i, 
                             treat.vec = treat.vec.in, give.treat=give.treat, treat.start=treat.start, nfix=nfix, N=N, do.clin.trial=do.clin.trial, start.trial = start.trial)
      
      #print(res.mf)
      
      all.mats.temp[, nfix + mf.c] <- res.mf
    }
    
    
    
    ####l1 and mf delay
    
    #take out last column to input into equations
    exp.delay.temp <- exposure.delay[, length(exposure.delay[1,])]
    mf.delay.temp <- mf.delay[, length(mf.delay[1,])]
    l1.delay.temp <- l1.delay
    
    exposure.delay[, inds.exp.mats] <- exposure.delay[, (inds.exp.mats -1)]
    mf.delay[, inds.mfd.mats] <- mf.delay[, (inds.mfd.mats - 1)] 
    
    #print(length(mf.delay.temp))
    #print(length(l1.delay.temp))
    
    
    #new L1 L2 and L3
    
    
    mf.temp <- rowSums(all.mats.cur[, (nfix+1) : (nfix + num.mf.comps)]) #sum mf over compartments, mf start in column 7
    
    
    all.mats.temp[, L1indi] <- calc.L1(beta, mf = mf.temp, mf.delay.in = mf.delay.temp, expos = tot.ex.ai, delta.vo=delta.vo,
                                       c.v=c.v, nuone=nuone, mu.v=mu.v, a.v=a.v, expos.delay = exp.delay.temp)
    all.mats.temp[, L2indi] <- calc.L2(nuone, L1.in = l1.delay.temp, mu.v, nutwo, mf = mf.delay.temp, a.v)
    all.mats.temp[, L3indi] <- calc.L3(nutwo, L2.in = all.mats.cur[, L2indi], a.H, g, mu.v, sigma.L0)
    
    #add new parasites 
    
    l1.delay <- all.mats.temp[, L1indi]
    mf.delay[, 1] <- rowSums(all.mats.cur[, (nfix+1) : (nfix + num.mf.comps)])
    exposure.delay[, 1] <- tot.ex.ai
    
    
    #new individual exposure for newborns, clear rows for new borns
    
    if(length(to.die) > 0)
    {
      #ex.vec[to.die] <- rlnorm(length(to.die), mu.ln, sigma.ln)
      ex.vec[to.die] <- rgamma(length(to.die), gam.dis, gam.dis)
      
      l.extras[to.die, ] <- 0 #establishing adult worms 
      mf.delay[to.die, 1] <- 0 #individuals dies so no contribution to l1s at this timestep
      l1.delay[to.die] <- 0
      
      treat.vec[to.die] <- NA
      
      all.mats.temp[to.die, cols.to.zero] <- 0 #set age, sex and parasites to 0 (includes L1, but not L2 L3)
      all.mats.temp[to.die, sexindi] <- rbinom(length(to.die), 1, 0.5) #draw sex
    }
    
    #update mf per snip and CMFL
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats.temp, nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    
    mean.mf.per.snip[i + 1] <- temp[[1]]
    
    CMFL[i + 1] <- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats.temp, ss.wt)
    
    prev.by.ss[i + 1] <- prevalence.for.age(age = age.prev.in, ss.in = temp, main.dat = all.mats.temp)
    
    #####################################################################################
    ## CTS-specific code - update clinical trial output vectors & apply loss to follow up
    #####################################################################################
    if (do.clin.trial & i>=(start.trial-1) & i<=stop.trial) { # changed to -1
      ## extract indicators for indiciduals still alive and in the cohort
      c1 <- which(all.mats.cur[,cohortindi]==1)
      c2 <- which(all.mats.cur[,cohortindi]==2)
      prev.by.ss.cohort[i+1,] <- c(length(which(temp[[2]][c1] > 0)) / length(c1),
                    length(which(temp[[2]][c2] > 0)) / length(c2))
      prev.diff[i+1] <- prev.by.ss.cohort[i+1, 1] - prev.by.ss.cohort[i+1, 2]
      ## infected or not output
      tmp <- as.numeric(temp[[2]]>0)
     # prev.by.ss.cohort[i+1,] <- c( mean(tmp[c1]), mean(tmp[c2]) )  
      sd.prev.by.ss.cohort[i+1,] <- c( sd(tmp[c1]), sd(tmp[c2]) )  
      mn.mf.by.ss.cohort[i+1,] <- c( mean(temp[[2]][c1]), mean(temp[[2]][c2]) )
      mf.diff[i+1] <- mn.mf.by.ss.cohort[i+1, 1] - mn.mf.by.ss.cohort[i+1, 2]
      sd.mf.by.ss.cohort[i+1,] <- c( sd(temp[[2]][c1]), sd(temp[[2]][c2]) )
   
      ## pull out a random number (if any) of c1 and c2 because of loss
      ## to follow up 
      ## number dropping out of cohorts
      cur.n <- sum(n.per.cohort[i,])
      n.drop <- rbinom(1,cur.n, prob.drop)
      ## which drop out
      if (n.drop>0) {
        drop <- sample.int(c(which(all.mats.cur[,cohortindi]==1), 
                         which(all.mats.cur[,cohortindi]==2)), n.drop)
        # update matrix ready for next time step
        all.mats.temp[drop,cohortindi] <- 0
      }
      
      n.per.cohort[i+1, ] <- c(sum(all.mats.temp[,cohortindi]==1), 
                               sum(all.mats.temp[,cohortindi]==2)) ## changed back to temp
      
    }
    #####################################################################################
    ## end of CTS chunk
    #####################################################################################
    
    all.mats[[i + 1]] <- all.mats.temp #update list of matrices
    
    L2[i + 1] <- mean(all.mats.temp[, L2indi])
    
    L1[i + 1] <- mean(all.mats.temp[, L1indi])
    
    #####parameters for fly larvae dynamics
  
    #debugging
    #l2par[i] <- mean(mf.delay.temp)
    #l2par2[i] <- mean(l1.delay.temp)
    
    ######
    ###bit all ove the place and in development!
    #####
    if(type.sero == 1)
      
    {  
      
      #age groups #automate this and have equal age group sizes?
      
      ###males 
      
      ind <- which(all.mats.temp[,ageindi] >= 0 & all.mats.temp[,ageindi] <= 2 & all.mats.temp[,sexindi] == 1)
      ind2 <- which(all.mats.temp[,ageindi] > 2 & all.mats.temp[,ageindi] <= 4 & all.mats.temp[,sexindi] == 1)
      ind3 <- which(all.mats.temp[,ageindi] > 4 & all.mats.temp[,ageindi] <= 6 & all.mats.temp[,sexindi] == 1)
      ind4 <- which(all.mats.temp[,ageindi] > 6 & all.mats.temp[,ageindi] <= 8 & all.mats.temp[,sexindi] == 1)
      ind5 <- which(all.mats.temp[,ageindi] > 8 & all.mats.temp[,ageindi] <= 10 & all.mats.temp[,sexindi] == 1)
      ind6 <- which(all.mats.temp[,ageindi] > 10 & all.mats.temp[,ageindi] <= 20 & all.mats.temp[,sexindi] == 1)
      ind7 <- which(all.mats.temp[,ageindi] > 20 & all.mats.temp[,ageindi] <= 30 & all.mats.temp[,sexindi] == 1)
      ind8 <- which(all.mats.temp[,ageindi] > 30 & all.mats.temp[,ageindi] <= 40 & all.mats.temp[,sexindi] == 1)
      ind9 <- which(all.mats.temp[,ageindi] > 40 & all.mats.temp[,ageindi] <= 50 & all.mats.temp[,sexindi] == 1)
      ind10 <- which(all.mats.temp[,ageindi] > 50 & all.mats.temp[,ageindi] <= 60 & all.mats.temp[,sexindi] == 1)
      ind11 <- which(all.mats.temp[,ageindi] > 60 & all.mats.temp[,ageindi] <= 70 & all.mats.temp[,sexindi] == 1)
      ind12 <- which(all.mats.temp[,ageindi] > 70 & all.mats.temp[,sexindi] == 1)
      
      inds.sero.males <- list(ind, ind2, ind3, ind4, ind5, ind6, ind7, ind8, ind9, ind10, ind11, ind12)  
      
      ###females 
      
      ind <- which(all.mats.temp[,ageindi] >= 0 & all.mats.temp[,ageindi] <= 2 & all.mats.temp[,sexindi] == 1)
      ind2 <- which(all.mats.temp[,ageindi] > 2 & all.mats.temp[,ageindi] <= 4 & all.mats.temp[,sexindi] == 1)
      ind3 <- which(all.mats.temp[,ageindi] > 4 & all.mats.temp[,ageindi] <= 6 & all.mats.temp[,sexindi] == 1)
      ind4 <- which(all.mats.temp[,ageindi] > 6 & all.mats.temp[,ageindi] <= 8 & all.mats.temp[,sexindi] == 1)
      ind5 <- which(all.mats.temp[,ageindi] > 8 & all.mats.temp[,ageindi] <= 10 & all.mats.temp[,sexindi] == 1)
      ind6 <- which(all.mats.temp[,ageindi] > 10 & all.mats.temp[,ageindi] <= 20 & all.mats.temp[,sexindi] == 1)
      ind7 <- which(all.mats.temp[,ageindi] > 20 & all.mats.temp[,ageindi] <= 30 & all.mats.temp[,sexindi] == 1)
      ind8 <- which(all.mats.temp[,ageindi] > 30 & all.mats.temp[,ageindi] <= 40 & all.mats.temp[,sexindi] == 1)
      ind9 <- which(all.mats.temp[,ageindi] > 40 & all.mats.temp[,ageindi] <= 50 & all.mats.temp[,sexindi] == 1)
      ind10 <- which(all.mats.temp[,ageindi] > 50 & all.mats.temp[,ageindi] <= 60 & all.mats.temp[,sexindi] == 1)
      ind11 <- which(all.mats.temp[,ageindi] > 60 & all.mats.temp[,ageindi] <= 70 & all.mats.temp[,sexindi] == 1)
      ind12 <- which(all.mats.temp[,ageindi] > 70 & all.mats.temp[,sexindi] == 0)
      
      inds.sero.females <- list(ind, ind2, ind3, ind4, ind5, ind6, ind7, ind8, ind9, ind10, ind11, ind12)  
      
      
      if(i == it.thresh) 
      {
        
        ######
        #####
        #initialize sero recording, depends on whether equilibrium is input or calculated during current run
        
        if(input.var == 0)
        {
          sero.prev.vec <- rep(0, N)
          ind <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sensitivity == 1) 
          sero.prev.vec[ind] <- 1
        }
        
        else #CHECK THIS
        {
          ind <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sero.prev.vec == 0 & sensitivity == 1) #which inds have worms but didn't at the last time step
          sero.prev.vec[ind] <- 1 #new (adult worm) infections become seropositive
          sero.prev.vec[to.die] <- 0 #newborns are negative
        }
        
        #########
        ##########
        
        ###############
        #age profiles for males and females
        ###############
        
        for(h in 1:length(inds.sero.males)) #could be vectorized
        {
          sprev.by.group.m[[h]] <- sum(sero.prev.vec[inds.sero.males[[h]]]) / length(inds.sero.males[[h]]) #males
          
          sprev.by.group.f[[h]] <- sum(sero.prev.vec[inds.sero.females[[h]]]) / length(inds.sero.females[[h]]) #females
          
          #mf prev by age 
          mfprev.by.group.m[[h]] <- length(which(temp[[2]][inds.sero.males[[h]]] > 0)) / length(inds.sero.males[[h]]) #males
          
          mfprev.by.group.f[[h]] <- length(which(temp[[2]][inds.sero.females[[h]]] > 0)) / length(inds.sero.females[[h]]) #females
        }
        
      }
      
      if(i > it.thresh)
      {
        ind <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sero.prev.vec == 0 & sensitivity == 1) #which inds have worms but didn't at the last time step
        sero.prev.vec[ind] <- 1 #new (adult worm) infections become seropositive
        sero.prev.vec[to.die] <- 0 #newborns are negative
        
        for(h in 1:length(inds.sero.males))
        {
          
          if(length(inds.sero.males[[h]]) <= 1) 
            
          {
            sprev.by.group.m[[h]] <- c(sprev.by.group.m[[h]], NA) #if no indivudals in age group, give NA
            mfprev.by.group.m[[h]] <- c(mfprev.by.group.m[[h]], NA)
          } 
          
          else
          {
            
            sprev.by.group.m[[h]] <- c(sprev.by.group.m[[h]], sum(sero.prev.vec[inds.sero.males[[h]]]) / length(inds.sero.males[[h]])) #males
            
            mfprev.by.group.m[[h]] <- c(mfprev.by.group.m[[h]], length(which(temp[[2]][inds.sero.males[[h]]] > 0)) / length(inds.sero.males[[h]]))
            
          }
          
          if(length(inds.sero.females[[h]]) <= 1) 
            
          {
            sprev.by.group.f[[h]] <- c(sprev.by.group.f[[h]], NA)
            mfprev.by.group.f[[h]] <- c(mfprev.by.group.f[[h]], NA)
          } 
          
          else
          {
            
            sprev.by.group.f[[h]] <- c(sprev.by.group.f[[h]], sum(sero.prev.vec[inds.sero.females[[h]]]) / length(inds.sero.females[[h]])) #females
            
            mfprev.by.group.f[[h]] <- c(mfprev.by.group.f[[h]], length(which(temp[[2]][inds.sero.females[[h]]] > 0)) / length(inds.sero.females[[h]]))
            
          }
          
        }
        
      }
      
      
      #stores sero prevalence for whole population and over 10s
      
      sero.prev[i] <- sum(sero.prev.vec) / N
      sero.prev.10[i] <- sum(sero.prev.vec[which(all.mats.cur[,ageindi] < 10)]) / length(which(all.mats.cur[,ageindi] < 10))
      
    } #closes serology
    
  } 
  ###############################################
  ## end of time loop
  ###############################################
  
  final.mat <- all.mats[[i + 1]] #save matrix from final time step before replacing with means
  
  #same process as at beginning of main for loop
  
  means <- vector(length=2)
  
  L3 <- sum(as.numeric(all.mats[[i]][, L3indi])) / N
  
  if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i]][, (nfix+1)])))} 
  else means[1] <- mean(as.numeric(rowSums(all.mats[[i]][, (nfix+1):(worms.start-1)]))) #sum up mf over compartments for each individual
  
  means[2] <- mean(as.numeric(rowSums(all.mats[[i]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )])))
  
  #############################################################################
  ## CTS-specific code - caclulate mean number of adult worms
  #############################################################################
  if (do.clin.trial) {
    means.cohort <- vector(length=2)
    means.cohort[1] <- means.cohort[2] <- 0
    if (i>=(start.trial-1) & i<=stop.trial) {  # changed to -1
      indi.tmp <- which(all.mats[[i]][,cohortindi]==1)
      if (length(indi.tmp)>0) {
      means.cohort[1] <- mean(as.numeric(rowSums(all.mats[[i]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ), drop=F])))
      } else {means.cohort[1]<-0}
      indi.tmp <- which(all.mats[[i]][,cohortindi]==2)
      if (length(indi.tmp)>0) {
      means.cohort[2] <- mean(as.numeric(rowSums(all.mats[[i]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ), drop=F])))
      } else {means.cohort[2]<-0}
    }
    all.mats[[i]] <- c(means, L3, means.cohort) #replace matrix at previous timestep with means for worms and L3
  } else {
    all.mats[[i]] <- c(means, L3) #replace matrix at previous timestep with means for worms and L3
  }
  ##############################################################################
  ## end of CTS-specific code chunck
  ##############################################################################
  
  means <- vector(length=2)
  
  L3 <- sum(as.numeric(all.mats[[i + 1]][, L3indi])) / N
  
  if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i + 1]][, (nfix+1)])))}
  else means[1] <- mean(as.numeric(rowSums(all.mats[[i + 1]][, (nfix+1):(worms.start-1)])))
  
  means[2] <- mean(as.numeric(rowSums(all.mats[[i + 1]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )])))
  
  #############################################################################
  ## CTS-specific code - caclulate mean number of adult worms
  #############################################################################
  if (do.clin.trial) {
    means.cohort <- vector(length=2)
    means.cohort[1] <- means.cohort[2] <- 0
    if (i>=(start.trial-1) & i<=stop.trial) { # changeed to -1
    indi.tmp <- which(all.mats[[i+1]][,cohortindi]==1)
    if (length(indi.tmp)>0){
      means.cohort[1] <- mean(as.numeric(rowSums(all.mats[[i+1]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ), drop=F])))
    } else {means.cohort[1]<-0} ## changed to fix for when there are no eligibles
    indi.tmp <- which(all.mats[[i+1]][,cohortindi]==2)
    if (length(indi.tmp)>0){
    means.cohort[2] <- mean(as.numeric(rowSums(all.mats[[i+1]][indi.tmp, worms.start : ((worms.start  + num.comps.worm * 3) - 1 ), drop=F])))
    } else {means.cohort[2]<-0}
    }
    all.mats[[i+1]] <- c(means, L3, means.cohort) #replace matrix at previous timestep with means for worms and L3
  } else {
    all.mats[[i+1]] <- c(means, L3) #replace matrix at previous timestep with means for worms and L3
  }
  ##############################################################################
  ## end of CTS-specific code chunck
  #####

  #ouput
  out.mf <- unlist(lapply(all.mats, FUN = `[`, 1))
  out.adults <- unlist(lapply(all.mats, FUN = `[`, 2))
  out.L3 <- unlist(lapply(all.mats, FUN = `[`, 3))
  
  #############################################################################
  ## CTS-specific code - mean number of adults output
  #############################################################################
  if (do.clin.trial) {
  out.adults.cohort.1 <- unlist(lapply(all.mats, FUN = `[`, 4))
  out.adults.cohort.2 <- unlist(lapply(all.mats, FUN = `[`, 5))
  }
  
  
  its <- seq(0:time.its)
  timez <- its * DT
  # which rows to output
  
  rows.out <- its[seq(1, length(its), by = round(max(every,DT)/DT))]
  
  ## would be good to add column labels to final data matrix
  if (get.equib==0) {
  out <- list(id.dat =final.mat, t = timez[rows.out], mn.mf = out.mf[rows.out], 
              mn.adult=out.adults[rows.out],mn.L3 = out.L3[rows.out], 
              mn.mf.ss = mean.mf.per.snip[rows.out], CMFL=CMFL[rows.out],
              pr.mf.ss =prev.by.ss[rows.out], sero.pr = sero.prev[rows.out],
              sero.pr.10 = sero.prev.10[rows.out], sero.pr.age.m = sprev.by.group.m[rows.out],
              pr.mf.age.m = mfprev.by.group.m[rows.out], sero.pr.age.f = sprev.by.group.f[rows.out],
              pr.mf.age.f = mfprev.by.group.f[rows.out], mn.L2=L2[rows.out], 
              mn.L1=L1[rows.out]) 
  ######################################################################
  ## CTS-specific code - modified output list
  ######################################################################
  if (do.clin.trial) {
    out <- c(out, list(applied.ss = applied.ss, 
                       prev.ss.cohort.1 = prev.by.ss.cohort[rows.out,1],
                       prev.ss.cohort.2 = prev.by.ss.cohort[rows.out,2],
                       mn.mf.ss.cohort.1 = mn.mf.by.ss.cohort[rows.out,1], 
                       mn.mf.ss.cohort.2 = mn.mf.by.ss.cohort[rows.out,2], 
                       sd.prev.by.cohort.1 = sd.prev.by.ss.cohort[rows.out,1],
                       sd.prev.by.cohort.2 = sd.prev.by.ss.cohort[rows.out,2],
                       sd.mf.ss.cohort.1 = sd.mf.by.ss.cohort[rows.out,1],
                       sd.mf.ss.cohort.2 = sd.mf.by.ss.cohort[rows.out,2],
                       n.cohort.1=n.per.cohort[rows.out,1],
                       n.cohort.2=n.per.cohort[rows.out,2],
                       mn.adult.cohort.1 = out.adults.cohort.1[rows.out], 
                       mn.adult.cohort.2 = out.adults.cohort.2[rows.out], 
                       mf.diff = mf.diff[rows.out], prev.diff = prev.diff[rows.out]))
  }
  ######################################################################
  ## end of CTS code chunk
  ######################################################################
  } else if (get.equib==1) {
    out <- list(final.mat = final.mat, ex.vec=ex.vec, sensitivity=sensitivity, sero.prev.vec=sero.prev.vec, 
                l.extras = l.extras, l1.delay=l1.delay, mf.delay=mf.delay, treat.vec.in=treat.vec.in, exposure.delay=exposure.delay)
  }
  out
  }

#save.image(paste("rfils/","eoibmout.",iter,".RData", sep=""))


