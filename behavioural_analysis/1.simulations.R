## simulate choices for the directed action task

library(R.matlab)
library(truncnorm)  # truncated normal
library(dplyr)
library(caroline)
library(tidyverse)

# Do you want to simulate online or offline data? (online meaning in scanner)
online = 1

# model variables
# first ones that change between versions
if (online == 1){
  NTr = 500 # no. of trials
} else if (online == 0){
  NTr = 1024
}
# now the consistent ones
set.seed(999)  # to set a seed value
N=21    # no. of ppts to simulate
Tsubj = rep(NTr,N) #vector of no. of trials for each ppt
Nopt = 2 # no. of options
Nrwd = 2 # no. of reward types
Qinits = matrix(0, nrow=2,ncol=2)
Vinits = rep(0.0,2)
Minits = rep(1,2)

# task variables
# first the ones that change
if (online == 1){
  tBack = 8 # no. of trials to look back for when assessing preference threshold
  rwdConsistent = 0.85 # token to dice is x% consistent
  states = 21
} else if (online == 0){
  tBack = 10 # no. of trials to look back for when assessing preference threshold
  rwdConsistent = 1.0 # token to dice is x% consistent
  states = 22
}
# now the consistent ones
H = 0.85
L = 0.15
threshold = 0.69
purseFull = 5

# task structure is the same but if online == 1 dont use first state
AR = c(H,L,H,H,0,H,H,L,L,0,H,H,H,0,H,0,0,0,0,0,0,0)
AG = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,H,H,H,H,0,H,H)
BR = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,H,H,L,L,0,H,H)
BG = c(L,H,L,L,0,H,H,H,H,0,H,H,L,0,H,0,0,0,0,0,0,0)
untilFull = c(0,0,0,1,0,0,0,0,2,0,0,0,1,0,0,0,0,0,2,0,0,0)
untilPrefer = c(1,2,1,0,0,2,0,2,0,0,1,0,0,0,2,1,0,1,0,0,2,0)
untilCount = c(20,20,20,0,8,20,4,20,0,8,20,4,0,8,20,20,4,20,0,8,20,4)
bankAfter = c(0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,1,0,0,0,0,2,0)
bankFull = c(3,3,3,2,2,2,3,3,1,1,1,3,2,2,2,2,3,3,1,1,1,3)

if (online == 1){
  AR = AR[2:22]
  AG = AG[2:22]
  BR = BR[2:22]
  BG = BG[2:22]
  untilFull = untilFull[2:22]
  untilPrefer = untilPrefer[2:22]
  untilCount = untilCount[2:22]
  bankAfter = bankAfter[2:22]
  bankFull = bankFull[2:22]
}

# variables to fill
rwd <- matrix(0,nrow=N,ncol=NTr)
rwd_devalued <- matrix(0,nrow=N,ncol=NTr)
choice <- matrix(0,nrow=N,ncol=NTr)
winCurr <- matrix(0,nrow=N,ncol=NTr)
isPurseFull_1 <- matrix(0,nrow=N,ncol=NTr)
isPurseFull_2 <- matrix(0,nrow=N,ncol=NTr)
purse = rep(0,2)
stateT <- matrix(0,nrow=N,ncol=NTr)
stateTimeT <- matrix(0,nrow=N,ncol=NTr)
purse1T <- matrix(0,nrow=N,ncol=NTr)
purse2T <- matrix(0,nrow=N,ncol=NTr)
ART <- matrix(0,nrow=N,ncol=NTr)
AGT <- matrix(0,nrow=N,ncol=NTr)
BRT <- matrix(0,nrow=N,ncol=NTr)
BGT <- matrix(0,nrow=N,ncol=NTr)
prefT <- matrix(0,nrow=N,ncol=NTr)
MR <- matrix(0,nrow=N,ncol=NTr)
MG <- matrix(0,nrow=N,ncol=NTr)
VA <- matrix(0,nrow=N,ncol=NTr)
VB <- matrix(0,nrow=N,ncol=NTr)
chosenNet <- matrix(0,nrow=N,ncol=NTr)
chosenDevalued <- matrix(0,nrow=N,ncol=NTr)
chosenRed <- matrix(0,nrow=N,ncol=NTr)
chosenGreen <- matrix(0,nrow=N,ncol=NTr)

### parameters
params <- read.csv('modelling/summaries/m10_summary.csv')
lr <- params$mean[7:28]
inv_temp <- params$mean[73:94]
M <- params$mean[29:50]

for (i in 1:N){
  
  v = Vinits
  Q = Qinits
  q = Vinits
  qc = Vinits
  m = Minits
  pe = Vinits
  pec = Vinits
  
  trialNr = 0
  state = 1
  stateTime = 0
  bank = 0
  purse = rep(0,2)
  
  for (t in 1:Tsubj[i]){
    
    ### rules governing state changes
    # which ever criteria is met first the state should reset
    trialNr = trialNr + 1
    stateTime = stateTime + 1

    
    # until count
    if (untilCount[state] > 0 && stateTime > untilCount[state]){
      stateTime = 1
      state = state + 1
      if (state > states){
        state = 1
        stateTime = 1
      }
    } 
    
    # until full
    if (untilFull[state] > 0 && purse[untilFull[state]] == purseFull){
      stateTime = 1
      state = state + 1
      if (state > states){
        state = 1
        stateTime = 1
      }
    }

    # until prefer
    # reset preference to 0
    # pref = 0
   
    # check if we're in a state that checks preference and we're above above no. of trials required to check pref
    if (untilPrefer[state] > 0 && stateTime > tBack){
      # check if we're in reinstatement states, here we want at least 12 trials
      if (state == 5 | state == 10 | state == 14 | state == 15 | state == 20){
        # we only want to check preference after 12 trials here
        if (stateTime > 11){
          pref = mean(choice[i,(t-1):(t-tBack)])-1
        } else {
          # manually set preferences before stateTime 12 so that it doesn't move to next state
          if (untilPrefer[state] == 1){
            pref = 1
          } else if (untilPrefer[state] == 2)
            pref = 0
        }
      } else {
        pref = mean(choice[i,(t-1):(t-tBack)])-1
      }
      # now check whether the preference needs to be inversed due to which dice we're checking for
      if (untilPrefer[state] == 1){
        pref = 1-pref
      } else if (untilPrefer[state] == 2){
        pref = pref
      }
      # now check if the preference is above threshold
      if (pref > threshold){
        stateTime = 1            # reset statetime to 1
        state = state + 1        # move to next state
        prefT[i,t] = pref        # record preference for this trial
        # check if moving to next state takes us over max state
        if (state > states){
          state = 1
          stateTime = 1
        }
      }

    # check if we're in a state that checks preference and we're above the no. of trials required before checking pref
    #if (untilPrefer[state] > 0 && stateTime > tBack){
    #  # check if were in reinstatement states, then we want to check after 12 trials
    #  if (state == 5 | state == 10 | state == 14
    #      | state == 15 | state == 20){
    #    if(stateTime < 12)
    #      pref = 999 # set pref to 999 if were in reinstatment and less than 12 trials
    #    else if (stateTime > 11){
    #      pref = 1-(mean(choice[i,(t-1):(t-tBack)])-1)
    #    }
    #  }
      # if it pref==999 then we want to skip preference check
    #  if (pref == 999){
    #    next
      # if not then check preference   
    #  } else if (pref != 999) {
    #    if (untilPrefer[state] == 1) {
    #      pref = 1-(mean(choice[i,(t-1):(t-tBack)])-1)
    #      if (pref > threshold){
    #        stateTime = 1
    #        state = state + 1
    #        prefT[i,t] = pref
    #        if (state == 22){
    #          state = 1
    #          stateTime = 1
    #        }
    #      }
    #    } else if (untilPrefer[state] == 2){
    #      pref = mean(choice[i,(t-1):(t-tBack)])-1
    #      if (pref > threshold){
    #        stateTime = 1
    #        state = state + 1
    #        prefT[i,t] = pref
    #        if (state == 22){
    #          state = 1
    #          stateTime = 1
    #        }
    #      }
    #    }
    #  }
    #}
    
    stateT[i,t] = state
    stateTimeT[i,t] = stateTime

    
    
    ## move purse to bank with the state change?
    # check every time we move to a new state
    if (state > 1 && stateTime == 1){
      if (bankAfter[state-1] == 1){
        purse[1] = 0
        bank = bank + 5
      }
      if (bankAfter[state-1] == 2){
        purse[2] = 0
        bank = bank + 5
      }
    }
    
    # reward probabilities
    rwdProbA = AR[state]+AG[state]
    rwdProbB = BR[state]+BG[state]
    
    ## Decision phase
    # set the motivation vector depending on whether the purse is full
    if (t > 1 && isPurseFull_1[i,t-1] == 1){
      m[1] = 1-M[i];
      m[2] = 1+M[i];
    } else if (t > 1 && isPurseFull_2[i,t-1] == 1){
      m[1] = 1+M[i];
      m[2] = 1-M[i];
    } else {
      m = Minits;
    }
    
    MR[i,t] <- m[1]
    MG[i,t] <- m[2]
     
    # collapse vectors and update v
    v[1] = (m[1]*Q[1,1])+(m[2]*Q[2,1]);
    v[2] = (m[1]*Q[1,2])+(m[2]*Q[2,2]);
    
    VA[i,t] <- v[1]
    VB[i,t] <- v[2]
    
    # generate choices
    choiceProb = 1/(1+exp(-inv_temp[i] * (v[1]-v[2])));
    choiceThreshold = runif(1)
    if (choiceProb > choiceThreshold){
      choice[i,t] = 1
    } else if (choiceProb < choiceThreshold){
      choice[i,t] = 2
    }
    
    # generate reward and token colour
    rewardThreshold = runif(1)
    # update reward for choice A
    if (choice[i,t] == 1){
      if (rewardThreshold < rwdProbA){
        rwd[i,t] = 1
        colourThreshold = runif(1)
        if (AR[state]>AG[state]){
          if (colourThreshold < rwdConsistent){
            winCurr[i,t] = 1
          } else if (colourThreshold > rwdConsistent){
            winCurr[i,t] = 2
          }
        } else if (AR[state]<AG[state]) {
          if (colourThreshold < rwdConsistent){
            winCurr[i,t] = 2
          } else if (colourThreshold > rwdConsistent){
            winCurr[i,t] = 1
          }
        } else if (rewardThreshold > rwdProbA) {
          rwd[i,t] = 0
      }
      }
    }
    
    
        
    # update reward for choice B
    if (choice[i,t] == 2){
      if (rewardThreshold < rwdProbB){
        rwd[i,t] = 1
        colourThreshold = runif(1)
        if (BG[state] > BR[state]){
          if (colourThreshold < rwdConsistent){
            winCurr[i,t] = 2
          } else if (colourThreshold > rwdConsistent){
            winCurr[i,t] = 1
          }
        } else if (BG[state] < BR[state]){
          if (colourThreshold < rwdConsistent){
            winCurr[i,t] = 1
          } else if (colourThreshold > rwdConsistent){
            winCurr[i,t] = 2
          }
        } else if (rewardThreshold > rwdProbB) {
          rwd[i,t] = 0
      }
      }
    }
    
    ## learning phase
    # update the value of the chosen option when a reward is received 
    # m during learning set to 1
    # counterfactual updating of unchosen
    
    
    if (rwd[i,t] == 1){
      pe[choice[i,t]] = rwd[i,t] - Q[winCurr[i,t],choice[i,t]];
      pec[3-choice[i,t]] = 0.0 - Q[winCurr[i,t],3-choice[i,t]];
      pec[choice[i,t]] = 0.0 - Q[3-winCurr[i,t],choice[i,t]];
      
      q[choice[i,t]] = 1 * (lr[i] * pe[choice[i,t]]);
      qc[3-choice[i,t]] = 1 * (lr[i] * pec[3-choice[i,t]]);
      qc[choice[i,t]] = 1 * (lr[i] * pec[choice[i,t]]);
      
      Q[winCurr[i,t],choice[i,t]] = Q[winCurr[i,t],choice[i,t]] + q[choice[i,t]];
      Q[winCurr[i,t],3-choice[i,t]] = Q[winCurr[i,t],3-choice[i,t]] + qc[3-choice[i,t]];
      Q[3-winCurr[i,t],choice[i,t]] = Q[3-winCurr[i,t],choice[i,t]] + qc[choice[i,t]];
      
    }
    
    # update value of the chosen option for both reward types when nothing is received 
    if (rwd[i,t] == 0){
      pe = rwd[i,t] - Q[,choice[i,t]];
      q = lr[i] * pe;
      Q[,choice[i,t]] = Q[,choice[i,t]] + q;
    }
    
    # bound Q values between 0 and 1
    for (j in 1:2){
      for (k in 1:2){
        if (Q[j,k] > 1){
          Q[j,k] = 1
          browser()
      } else if (Q[j,k] < 0){
        Q[j,k] = 0
        browser()
      }
      }
    }
    
    
    ART[i,t] = Q[1,1]
    AGT[i,t] = Q[2,1]
    BRT[i,t] = Q[1,2]
    BGT[i,t] = Q[2,2]
    
    chosenNet[i,t] = Q[1,choice[i,t]] + Q[2,choice[i,t]];
    chosenDevalued[i,t] = (m[1]*Q[1,choice[i,t]]) + (m[2]*Q[2,choice[i,t]]);
    chosenRed[i,t] = Q[1,choice[i,t]]
    chosenGreen[i,t] = Q[2,choice[i,t]]
    
    ## update purses
    if (winCurr[i,t] == 1){
      purse[1] = purse[1] + 1
    } else if (winCurr[i,t] == 2)
      purse[2] = purse[2] + 1
    
    if (purse[1] == purseFull | purse[1] > purseFull){
      purse[1] = purseFull
      isPurseFull_1[i,t] = 1
    }
    
    if (purse[2] == purseFull | purse[2] > purseFull){
      purse[2] = purseFull
      isPurseFull_2[i,t] = 1
    }
    
    purse1T[i,t] = purse[1]
    purse2T[i,t] = purse[2]
    
    ## move purse to bank? 
    # bank when full? 
    if (bankFull[state] == 3 && rwd[i,t] == 1){
      if (purse[winCurr[i,t]] == purseFull){
        purse[winCurr[i,t]] = 0
        bank = bank + 5
        isPurseFull_1[i,t] = 0
        isPurseFull_2[i,t] = 0
      }
    }
    if (bankFull[state] == 2){
      if (purse[2] == purseFull){
        purse[2] = 0
        bank = bank + 5
        isPurseFull_2[i,t] = 0
      }
    }
    if (bankFull[state] == 1){
      if (purse[1] == purseFull){
        purse[1] = 0
        bank = bank + 5
        isPurseFull_1[i,t] = 0
      }
    }
    
    if (winCurr[i,t] > 0){
      rwd_devalued[i,t] = rwd[i,t]*m[winCurr[i,t]]
    } else {
      rwd_devalued[i,t] = 0
    }
    
  }
  }
}



## organise data for behavioural analysis 
# 1 column for each variable 

data <- as.data.frame(choice) %>%
  pivot_longer(everything()) %>%
  rename(chosenDice = value)

x <- as.data.frame(rwd) %>%
  pivot_longer(everything()) %>%
  rename(totalWin = value)

data <- cbind(data,x)

x <- as.data.frame(stateT) %>%
  pivot_longer(everything()) %>%
  rename(state = value)

data <- cbind(data,x)

x <- as.data.frame(MR) %>%
  pivot_longer(everything()) %>%
  rename(MR = value)

data <- cbind(data,x)

x <- as.data.frame(MG) %>%
  pivot_longer(everything()) %>%
  rename(MG = value)

data <- cbind(data,x)

x <- as.data.frame(VA) %>%
  pivot_longer(everything()) %>%
  rename(VA = value)

data <- cbind(data,x)

x <- as.data.frame(VB) %>%
  pivot_longer(everything()) %>%
  rename(VB = value)

data <- cbind(data,x)

x <- as.data.frame(chosenNet) %>%
  pivot_longer(everything()) %>%
  rename(chosenNet = value)

data <- cbind(data,x)

x <- as.data.frame(chosenDevalued) %>%
  pivot_longer(everything()) %>%
  rename(chosenDevalued = value)

data <- cbind(data,x)

x <- as.data.frame(rwd_devalued) %>%
  pivot_longer(everything()) %>%
  rename(rwd_devalued = value)

data <- cbind(data,x)

x <- as.data.frame(winCurr) %>%
  pivot_longer(everything()) %>%
  rename(winCurr = value)

data <- cbind(data,x)

x <- as.data.frame(chosenRed) %>%
  pivot_longer(everything()) %>%
  rename(chosenRed = value)

data <- cbind(data,x)

x <- as.data.frame(chosenGreen) %>%
  pivot_longer(everything()) %>%
  rename(chosenGreen = value)

data <- cbind(data,x)

x <- as.data.frame(stateTimeT) %>%
  pivot_longer(everything()) %>%
  rename(stateTime = value)

data <- cbind(data,x)

x <- as.data.frame(isPurseFull_1) %>%
  pivot_longer(everything()) %>%
  rename(isPurseFull_1 = value)

data <- cbind(data,x)

x <- as.data.frame(isPurseFull_2) %>%
  pivot_longer(everything()) %>%
  rename(isPurseFull_2 = value)

data <- cbind(data,x)

x <- as.data.frame(prefT) %>%
  pivot_longer(everything()) %>%
  rename(prefT = value)

data <- cbind(data,x)

data <- data %>%
  select(chosenDice,totalWin,state,stateTime,MR,MG,isPurseFull_1,isPurseFull_2,VA,VB,chosenNet,chosenDevalued,rwd_devalued, winCurr, chosenRed, chosenGreen, prefT) %>%
  mutate(vDiff = chosenNet - chosenDevalued) 
data['subject_nr'] <- rep(1:21, each = NTr)
data['allTrialIndex'] <- rep(1:NTr, 21)

# save the data 
if (online==0){
  write.csv(data, 'raw_data/simulations/simulated_choice_data_FM_determ2.csv')
} else if (online==1){
  write.csv(data, 'raw_data/simulations/simulated_choice_data_M.csv')
}

cor(data$chosenNet,data$chosenDevalued,method='pearson')
cor(data$chosenNet[1:NTr],data$chosenDevalued[1:NTr],method='pearson')

Devaluedstate_corr <- data %>%
  filter(MR == 0 | MG == 0) 
cor(Devaluedstate_corr$chosenNet,Devaluedstate_corr$chosenDevalued,method='pearson')

cor(data$totalWin,data$rwd_devalued,method='pearson')

y <- data$chosenNet - data$chosenDevalued
trialNumber <- 1:NTr
diff <- as.data.frame(trialNumber)
diff['vDiff'] <- y[1:NTr]

# get the states and statetimes of high signal
vDiffTime <- data %>%
  select(state,stateTime,allTrialIndex,subject_nr,vDiff) %>%
  filter(vDiff > 0.001)

stateAVG <- vDiffTime %>%
  select(state,stateTime,vDiff) %>%
  group_by(state) %>%
  mutate(stateMean = mean(vDiff)) %>%
  distinct(state,stateMean)
  

vDiffLine <- ggplot(diff, aes(y=vDiff, x=trialNumber)) +
  geom_line(colour='red')

vDiffLine <- ggplot(data, aes(y=signal, x=trialNr)) +
  geom_line(colour='red') +
  theme_bw() +
  labs(title='Devaluation signal across trials in the two currencies task', x='Trial Number', y='Devaluation signal')

sig = mean(diff$vDiff > 0.001)

sig = mean(data$signal > 0.001)

RG_corr = cor(data$chosenGreen,data$chosenRed,method='pearson')

## plot correlations and signal levels across across degrees of uncertainty 

level_0.65 = c(RG_corr,sig)

write.csv(data, 'simulated_choice_data.csv')

### Turn chosenNet chosenDevalued difference signal into a simulated HRF
# get real trial times

times <- read.csv("jittered_ITI_6.csv") %>%
  select(startDice) %>%
  rowid_to_column()

times['time'] <- 0

for (i in times$rowid){
  if (times[i,]$rowid == 1) next
  times[i,]$time = times[i,]$startDice - times[i-1,]$startDice
}
  
times <- times %>%
  mutate(cumulativeTime = round(cumsum(time))) 

times['signal'] <- diff$vDiff

totalTime = times$cumulativeTime[NTr]

# remove all signal that isn't of meaningful magnitude 
for (i in times$rowid){
  if (times[i,]$signal < 0.001){
    times[i,]$signal = 0
  }
}

# make vector 1 row for each second
time <- 1:totalTime

# create the HRF data frame
hrf <- as.data.frame(time) %>%
  rowid_to_column()
hrf['trials'] <- 0
hrf['signal'] <- 0

# set trial to 1 when a trial occurred on that second 
for (i in hrf$time){
  for (j in times$rowid){
    if (hrf[i,]$time == times[j,]$cumulativeTime){
      hrf[i,]$trials = 1
    }
  }
}

# alternatively set trials to 1 when signal is present
for (i in hrf$time){
  if(hrf[i,]$signal > 0){
    hrf[i,]$trials = 1
  }
}

# set signal to corresponding time point
for (i in hrf$time){
  for (j in times$rowid){
    if (hrf[i,]$time == times[j,]$cumulativeTime){
      hrf[i,]$signal = times[j,]$signal
    }
  }
}

write.csv(hrf, 'hrf_simulation_100_80.csv')

signal_df <- as.data.frame(time)
signal = hrf$signal
signal_df['signal'] <- signal

write.csv(signal_df, 'signal_2.csv')

## actual value instead of EV

actual_time <- read.csv('jittered_ITI_6.csv') %>%
  select(startReward,startDice) %>%
  rowid_to_column() %>%
  mutate(diff=startReward-startDice)

actual_time['time'] <- times$cumulativeTime

actual_time <- actual_time %>%
  mutate(RewardTime=time+diff)

actual_signal <- data %>%
  select(totalWin, rwd_devalued) %>%
  mutate(actual_diff = totalWin - rwd_devalued)

actual_diff = actual_signal$actual_diff[1:NTr]

totalTime = actual_time$RewardTime[NTr]

time_2 = 1:totalTime

hrf_2 <- as.data.frame(time_2) %>%
  rowid_to_column()

hrf_2['trials'] <- 0
hrf_2['actual_signal'] <- 0
  
actual_time['actual_signal'] <- actual_diff

hrf['cumulativeTime_actaul'] <- actual_time$RewardTime

actual_time <- actual_time %>%
  mutate(time = round(RewardTime)) 

# set signal to corresponding time point
for (i in hrf_2$time_2){
  for (j in actual_time$rowid){
    if (hrf_2[i,]$time_2 == actual_time[j,]$time){
      hrf_2[i,]$actual_signal = actual_time[j,]$actual_signal
    }
  }
}

# set time to 1 when trial occurs
for (i in hrf_2$time_2){
  for (j in actual_time$rowid){
    if (hrf_2[i,]$time_2 == actual_time[j,]$cumulativeTime_actual){
      hrf_2[i,]$trials = 1
    }
  }
}

# check correlations

cor(hrf$signal,hrf$actual_signal,method='pearson')
cor(hrf$signal,hrf$PE,method='pearson')

### Now lets do the same for chosenGreen chosenRed correlation
# get real trial times

times['chosenRed'] <- data$chosenRed[1:NTr]
times['chosenGreen'] <- data$chosenGreen[1:NTr]

times <- times %>%
  mutate(chosen_diff = chosenRed-chosenGreen)

hrf['chosenRed'] <- 0
hrf['chosenGreen'] <- 0

# set chosenRed and Green to corresponding time point
for (i in hrf$time){
  for (j in times$rowid){
    if (hrf[i,]$time == times[j,]$cumulativeTime){
      hrf[i,]$chosenRed = times[j,]$chosenRed
      hrf[i,]$chosenGreen = times[j,]$chosenGreen
    }
  }
}
# how red is your choice
hrf <- hrf %>%
  mutate(Chosen_diff = chosenRed-chosenGreen)

write.csv(hrf, 'hrf_chosen_colour.csv')

signal_df <- as.data.frame(time)
signal = hrf$signal
signal_df['signal'] <- signal

write.csv(signal_df, 'signal_2.csv')

# add reward and net chosen EV obb to hrf df

hrf['reward'] <- 0
hrf['chosenEV'] <- 0

# set reward to reward time point
for (i in hrf$time){
  for (j in times$trialNr){
    if (hrf[i,]$time == times[j,]$rewardTime){
      hrf[i,]$reward=times[j,]$netReward
    }
  }
}

# set ev to decision time point
for (i in hrf$time){
  for (j in times$trialNr){
    if (hrf[i,]$time == times[j,]$cumulativeTime){
      hrf[i,]$chosenEV=times[j,]$netChosen
    }
  }
}

## other regressors 

# trialNr, just rename rowid here
times <- times %>%
  rename(trialNr=rowid)

# net chosen expected value, red+green for the chosen option

times <- times %>%
  mutate(netChosen = chosenRed+chosenGreen)

# net reward received

a <- data$totalWin[1:NTr]

times['netReward'] <- a

# RPE

times <- times %>%
  mutate(RPE = netReward - netChosen)

# append actual signal and reward time to times df

times['rewardTime'] <- actual_time$time
times['actual_signal'] <- actual_time$actual_signal

## adapt times df for making regressors netChosen, netReward, and RPE sub + obb

times['chosenDevalued'] <- data$chosenDevalued[1:NTr]
times['rwdDevalued'] <- data$rwd_devalued[1:NTr]

times <- times %>%
  mutate(chosenNetSO = netChosen + chosenDevalued) %>%
  mutate(rwdSO = netReward + rwdDevalued)



### make text files for efficiency analysis 

# ev signal

ev_signal_text <- times %>%
  select(cumulativeTime,signal)
ev_signal_text['zeros'] <- 0
ev_signal_text <- ev_signal_text %>%
  select(cumulativeTime,zeros,signal)

write.table(ev_signal_text, file='ev_signal_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# ov signal

ov_signal_text <- times %>%
  select(rewardTime,actual_signal)
ov_signal_text['zeros'] <- 0
ov_signal_text <- ov_signal_text %>%
  select(rewardTime,zeros,actual_signal)

write.table(ov_signal_text, file='ov_signal_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# chosen colour difference 

chosen_diff_text <- times %>%
  select(cumulativeTime,chosen_diff)
chosen_diff_text['zeros'] <- 0
chosen_diff_text <- chosen_diff_text %>%
  select(cumulativeTime,zeros,chosen_diff)

write.table(chosen_diff_text, file='chosen_diff_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# RPE Signal

RPE_text <- times %>%
  select(cumulativeTime,RPE)
RPE_text['zeros'] <- 0
RPE_text <- RPE_text %>%
  select(cumulativeTime,zeros,RPE)

write.table(RPE_text, file='RPE_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# trials 

trials_text <- times %>%
  select(cumulativeTime)
trials_text['zeros'] <- 0
trials_text['trials'] <- 1
trials_text <- trials_text %>%
  select(cumulativeTime,zeros,trials)

write.table(trials_text, file='trials_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# net chosen EV

chosenEV_text <- times %>%
  select(cumulativeTime,netChosen)
chosenEV_text['zeros'] <- 0
chosenEV_text <- chosenEV_text %>%
  select(cumulativeTime,zeros,netChosen)

write.table(chosenEV_text, file='chosenEV_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# net RWD

reward_text <- times %>%
  select(rewardTime,netReward)
reward_text['zeros'] <- 0
reward_text <- reward_text %>%
  select(rewardTime,zeros,netReward)

write.table(reward_text, file='reward_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# net chosen EV sub + ob

chosenEV_SO_text <- times %>%
  select(cumulativeTime,chosenNetSO)
chosenEV_SO_text['zeros'] <- 0
chosenEV_SO_text <- chosenEV_SO_text %>%
  select(cumulativeTime,zeros,chosenNetSO)

write.table(chosenEV_SO_text, file='chosenEV_SO_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# net RWD sub + ob

reward_SO_text <- times %>%
  select(rewardTime,rwdSO)
reward_SO_text['zeros'] <- 0
reward_SO_text <- reward_SO_text %>%
  select(rewardTime,zeros,rwdSO)

write.table(reward_SO_text, file='reward_SO_85_85_final.txt',sep="\t", row.names=FALSE,col.names=FALSE)

# save dfs as csv

write.csv(times, 'simulations_85_85.csv')
