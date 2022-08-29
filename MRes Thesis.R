# NOTATION
# s_0: prior belief about successes on treatment D before Phase II.
# f_0: prior belief about failures on treatment D before Phase II.
# s_1: number of successes after Phase II for treatment D.
# f_1: number of failures after Phase II for treatment D.
# s_1C: prior belief about successes on treatment C before Phase III
# f_1C: prior belief about failures on treatment C before Phase III
# p_C: proportion of successes for treatment C Phase II
# C: cost of not developing drug.
# c: cost of recruiting a patient.
# Pop: Population benefiting from treatments.
# n: patients in Phase III.
# R: Revenue from the new treatment (per patient).
# T_a: patients for the selected Phase II.
# theta_C: quantity indicating effectiveness of treatment C before Phase II.
# alpha: significance level for Phase II.
# beta: type II error for Phase II
# k: minimal detectable difference
# n_1: list of Phase II samples


library(hash)
library(Exact)
library(VGAM)
CRDP <- function(R,Pop,n_1,alpha,beta,theta_C,k, p_C, C,c, s_0, f_0, s_1C, f_1C){
  #Phase II possible outcomes
  #Creating a list of list to store all the possible outcomes from Phase II
  Outcomes1 = list(list()) 
  ind = 1
  for (n1 in n_1){
    Outcomes2 = list(list())
    i = 1
    #For every number of successes apo to n1 
    for (s_1 in 0:n1){
      #Failures are "deterministic" once you know n1, s_1
      f_1 = n1 - s_1
      #Append successes and failures to a temporary list
      Outcomes2[[i]] = list(c(s_1,f_1))
      i = i + 1
    }
    #Store the list for n1 to the specific index and go to the next one
    Outcomes1[[ind]] = Outcomes2
    ind = ind + 1
  }
  
  
  
  #use library(hash) to introduce dictionaries
  #Here we will create a list of dictionaries for every possible phase II
  #By the end list_of_phase_II will have as key (s1,f1) and as value the sample size n for Phase III
  list_of_phase_II = list()
  list_of_samples = hash()
  maximum_value = 0
  maximum_successes = 0
  for (i in 1:length(Outcomes1)){
    for (j in 1:length(Outcomes1[[i]])){
      #s1,f1
      successes = Outcomes1[[i]][[j]][[1]][1]
      failures = Outcomes1[[i]][[j]][[1]][2]
      
      #create dictionary with successes and failures concatenated as key, and the sample n as value
      key = paste(successes,",",failures)
      recovered_sample_size = successes + failures
      prop_successes = successes / (recovered_sample_size)

      #Find sample size for phase III for each possible outcome of Phase II
      
      n = floor((qnorm(1-alpha) + qnorm(1-beta))^2 *
                  (prop_successes * (1 - prop_successes) +p_C*(1-p_C))/k^2)
      
      
      #Both ifs are here for dimension purposes
      if (n > maximum_value){
        maximum_value = n
      }
      if (successes > maximum_successes){
        maximum_successes = successes
      }
      list_of_samples[[key]] = n
    }
    list_of_phase_II[[i]] = list_of_samples
    list_of_samples = hash()
  }

  
  #Reward has s_2D, s_2C, s1, f1, k(number of phase II trials)
  Reward2<- array(0, dim = c(floor(maximum_value/2) + 1,
                             floor(maximum_value/2) + 1,
                             maximum_successes+1,
                             maximum_successes+1,
                             length(list_of_phase_II)))
  
  
  
  #Here we populate Reward at t = 2 (only if we had go decision)
  #We start by taking every possible Phase II (k)
  for (k in 1:length(list_of_phase_II)){
    #For every (s1,f1) 
    for (key in ls(list_of_phase_II[[k]])){
      sum_go = 0
      sum_no_go = 0
      #Find s1
      s_1 = as.numeric(strsplit(key, ",")[[1]][1])
      #Find f1
      f_1 = as.numeric(strsplit(key, ",")[[1]][2])
      #Find n
      n = list_of_phase_II[[k]][[key]]
      #Find where we reject H0
      reject_region = exact.reject.region(floor(n/2),floor(n/2),
                                          "greater", alpha = alpha,method = "Fisher")
      #Now for every possible combination (given k,s1,f1)
      for (s_2D in 0:(floor(n/2))){
        for (s_2C in 0:(floor(n/2))){
          s_1D = s_0 + s_1
          f_1D = f_0 + f_1
          expected_prob_C = dbetabinom.ab(s_2C, n/2, s_1C,f_1C)
          expected_prob_D = dbetabinom.ab(s_2D ,n/2, s_1D,f_1D)
          #Find reward at t = 2, taken from R_2 MDP
          Reward2[s_2D + 1, s_2C + 1, s_1+1,f_1+1,k] = (s_2D + s_2C + (reject_region[s_2D+1,s_2C+1])*R*Pop-
                                                          (1- reject_region[s_2D+1,s_2C+1] * C * Pop))
        }
      }
    }
  }
  #We shall go back to t = 1 and populate Rewards for go and no go there
  Reward1_go = array(0, dim = c(dim(Reward2)[1],
                                dim(Reward2)[2],
                                dim(Reward2)[3],
                                dim(Reward2)[4],
                                dim(Reward2)[5]))
  Reward1_no_go = array(0, dim = c(dim(Reward2)[1],
                                   dim(Reward2)[2],
                                   dim(Reward2)[3],
                                   dim(Reward2)[4],
                                   dim(Reward2)[5]))
  
  #Similar to before, for every Phase II
  for (k in 1:length(list_of_phase_II)){
    #For every (s1,f1)
    for (key in ls(list_of_phase_II[[k]])){
      s_1 = as.numeric(strsplit(key, ",")[[1]][1])
      f_1 = as.numeric(strsplit(key, ",")[[1]][2])
      n = list_of_phase_II[[k]][[key]]
      
      #For every possible combination s2D, s2C (given s1,f1,k)
      for (s_2D in 0:(floor(n/2))){
        for (s_2C in 0:(floor(n/2))){
          s_1D = s_0 + s_1
          f_1D = f_0 + f_1
          #Find posterior probability of reaching that specific combination
          expected_prob_C = dbetabinom.ab(s_2C, n/2, s_1C,f_1C)
          expected_prob_D = dbetabinom.ab(s_2D ,n/2, s_1D,f_1D)
          
          #Find every reward at t=1 for every possible s_2D, s_2C given s1, f1, k
          if (Reward2[s_2D+1,s_2C+1,s_1+1,f_1+1,k] != 0){
            #We use product of the two as we assume independence
            Reward1_go[s_2D+1,s_2C+1,s_1+1,f_1+1,k] = expected_prob_C*expected_prob_D*(Reward2[s_2D+1,s_2C+1,s_1+1,f_1+1,k]) + s_1 - n*c
            Reward1_no_go[s_2D+1,s_2C+1,s_1+1,f_1+1,k] = s_1 - C*Pop
          }
          else{
            Reward1_go[s_2D+1,s_2C+1,s_1+1,f_1+1,k] = 0
            Reward1_no_go[s_2D+1,s_2C+1,s_1+1,f_1+1,k] = 0
          }
        }
      }
    }
  }#
  
  Action1 = list()
  #The output here will be a list of dictionaries
  #where the key will be (s1,f1,cumulative reward of optimal policy), for every possible s1,f1 given k
  #and the value will be either 1 (no go) or 2 (go)
  
  #For every possible Phase II outcome
  for (k in 1:length(list_of_phase_II)){
    #We create a dictionary for each k
    dictionary_of_decision = hash()
    for (key in ls(list_of_phase_II[[k]])){
      sum_go = 0
      sum_no_go = 0
      s_1 = as.numeric(strsplit(key, ",")[[1]][1])
      f_1 = as.numeric(strsplit(key, ",")[[1]][2])
      n = list_of_phase_II[[k]][[key]]
      sum_go = sum(Reward1_go[,,s_1+1,f_1+1,k])
      sum_no_go = sum(Reward1_no_go[,,s_1+1,f_1+1,k])
      
      #Choose go if weighted sum is larger than no go
      if((sum_go >= sum_no_go)&(sum_go > 0)){
        key = paste(s_1,",",f_1,",",sum_go)
        dictionary_of_decision[[key]] = 2
      }
      if((sum_go >= sum_no_go)&(sum_go<0)){
        sum_go = 0
        key = paste(s_1,",",f_1,",",sum_go)
        dictionary_of_decision[[key]] = 1
      }
      if(sum_go < sum_no_go){
        key = paste(s_1,",",f_1,",",sum_no_go)
        dictionary_of_decision[[key]] = 1
      }
    }
    #For phase k we have the dictionary
    Action1[[k]] = dictionary_of_decision
  }
  
  #Here we iterate over every phase II trial
  #The goal is to find which phase II trial returns the largest expected weighted reward
  
  values = c()
  #For every k
  for (k in 1:length(list_of_phase_II)){
    sum_t_0 = 0
    #For every (s1,f1,reward), given we chose trial k
    for (key_reward in ls(Action1[[k]])){
      s_1 = as.numeric(strsplit(key_reward, ",")[[1]][1])
      f_1 = as.numeric(strsplit(key_reward, ",")[[1]][2])
      reward = as.numeric(strsplit(key_reward, ",")[[1]][3])
      T_a = s_1 + f_1
      expected_prob = dbetabinom.ab(s_1, T_a, s_0, f_0)
      sum_t_0 = sum_t_0 + expected_prob*reward - T_a*(theta_C + c)
    }
    values = c(values,sum_t_0)
  }
  return(list(Action1 = Action1, values = values))
}

