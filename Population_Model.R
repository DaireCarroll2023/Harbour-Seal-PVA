#################################
#Daire Carroll, University of Gothenburg, 2021, daire.carroll@bioenv.gu.se
#################################
#function to calculate the age structured population at time t given starting population P
#Parameters:
#prarameters
#Ka = population carrying capacity...  numeric (vector for multiple populations)
#t = number of cycles (years)... numeric
#A = max age... numeric
#B = Birthrates, missing value will be 0...   vector
#S = Survival, , missing value will be 0...   vector
#P = initial population structure, can be established based on the dominant eigen value of the leslie matrix based on structure function and population size N...  vector for one population or matrix with each column representing a population

structure = function(B,S,N){
  L = matrix(B,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev = as.numeric(ev$vectors[,1])
  P = c()
  for(i in 1:A){
    P[i] = ev[i]/sum(ev)*N
  } 
  return(P)
} #using stable age structure to estimate number of seals in each age class

#H = Proportion of population lost due to hunting, can be 0... numeric
#H_delay = 0 Numnber of cycles before hunting is implemented, can be 0... numeric
#Target_ages = Ages targeted with disproportionate level of hunting, can be NULL... vector
#H_skew = The proportion of total hunting to be assigned to target age classes, remaining hunting will be distributed among other age classes based on proportional size. If one of the target classes is <1, hunting will be distributed among remaining Target_ages, if all are <1, it will be distributed among remaining age classes... numeric 
#H_d = An alternative method for skewing hunting, the proportion of total hunting to be allocated in each class... vector

#E = Entanglement, currently a single value of seals removed from the each population per cycle, will be modified to have optional age class skew or be a function of populations size... numeric

#f = Modification to birthrates, will reduce initial fertility down by f % until 0, f can be 0... numeric of a vector of additive effects 

#R = level of random variation per cycle (0.05 = random value within 5% +- calculated value based on beta distribution a = b = 5)... numeric

#P_D = Probability of a disease outbreak in each cycle... numeric
#O_D = Outcome of disease, maximum proportion of an age class killed during an outbreak... numeric
#W_D = Proportion of O_D experienced by age class... vector

#M = Migration: proportion of a population which is present in each population at the end of each cycle... matrix 

################################

PopMod = function(t,A,Ka,P,B,S,f,H,H_delay,E,R,Target_ages,H_skew,H_d,P_D,O_D,W_D,M){
  
  null_args = list(Target_ages,H_skew,H_d,P_D,O_D,W_D,M)
  null_args_names = c("Target_ages","H_skew","H_d","P_D","O_D","W_D","M")
  for(i in 1:length(null_args)){
    if(is.null(null_args[[i]])){
      assign(paste0(null_args_names[i]), 0)
    }
  }
  
  if(nrow(P)!=A || (length(S))!=A || length(B)!=A){
    stop('Lengths of P, S, & B must = A')
  }
  
  if(t==0){		
    return(P)
  }
  
  L = matrix(B,ncol = A) #begin to define the leslie matrix, to be populated
  
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }	#populate leslie matrix with survival values from S
  
  ev = eigen(L)
  ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
  
  f = sum(f)	#allow for multiple contributions to modifications of birthrate
  L[1,] = L[1,] - f*L[1,] #apply constant fertility modifier
  
  PGK = function(P){ 
    
    P_out = matrix(0, nrow = nrow(P),ncol = ncol(P))
    
    check = function(Pop){
      for(i in 1:length(Pop)){
        if(Pop[i] < 0){
          Pop[i] = 0
        }
      }	 
      return(Pop)
    } #remove any negative population values
    
    Keffect = function(K,ev1,TP){ 
      
      ((K+(ev1-1)*TP)/K)^-1
      
    } #calculate the effect of carrying capactity at t based on population size at t-1
    
    apply_H = function(PI){
      
      basic_H = function(PI,H_v){
        
        PI = PI-PI*H_v
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      even_spread = function(reduction,neg_v,pos_v){
        
        while(length(neg_v )>0){
          
          deficit = c()
          
          for(i in 1:length(neg_v)){
            deficit[i] = reduction[neg_v[i]]
          }
          
          deficit = sum(deficit)
          
          for(i in 1:length(neg_v)){
            
            reduction[neg_v[i]] = 0
            
          }
          
          reduction_L = sum(reduction)
          
          for(i in 1:length(pos_v)){
            
            reduction[pos_v[i]] = reduction[pos_v[i]] + deficit*(reduction[pos_v[i]]/reduction_L) 
            
            neg_v = which(reduction<0)
            pos_v = which(reduction>0)
            
          }
          
        }
        
        return(reduction)
        
      } #when uneven hunting leads to depletion of an age class, deficit will be shared out based on proportion of the population 
      
      if(H_skew == 0 && H_d == 0){
        
        PI = basic_H(PI,H)
        
        return(PI) #with no hunting skew, basic hunting is applied across classes
        
      }else if(H_skew == 0){
        
        reduction = PI-H_d*H*sum(PI)
        
        neg_v = which(reduction<0)
        pos_v = which(reduction>0)
        
        PI = even_spread(reduction,neg_v,pos_v)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }else if(H_skew > 0){
        
        reduction = rep(NA, times = A)
        
        for(i in 1:length(Target_ages)){
          
          reduction[Target_ages[i]] =  PI[Target_ages[i]] - H*H_skew*sum(PI)
          
        }
        
        for(i in 1:length(reduction)){
          if(is.na(reduction[i])==TRUE){
            reduction[i] = PI[i] - (1-H_skew)*H*(PI[i]/sum(PI)) 
          }			
        }
        
        neg_v = which(reduction<0)
        pos_v = which(reduction>0)
        
        while(length(neg_v) > 0 && length(neg_v) < length(Target_ages)){
          
          deficit = c()
          
          for(i in 1:length(neg_v)){
            deficit[i] = reduction[neg_v[i]]
          }
          
          deficit = sum(deficit)
          
          for(i in 1:length(neg_v)){
            
            reduction[neg_v[i]] = 0
            
          }
          
          spread = c()
          
          for(i in 1:length(Target_ages)){
            
            if(reduction[Target_ages[i]]>0){
              
              spread[i] = Target_ages[i]
              
            }
            
          }
          
          for(i in 1:length(spread)){
            
            reduction[spread[i]] = reduction[spread[i]] + deficit/length(spread)
            
          }
          
          neg_v = which(reduction<0)
          pos_v = which(reduction>0)
          
          PI = reduction
          
        } #the deficit is spread between the targeted age classes
        
        PI = even_spread(reduction,neg_v,pos_v)
        
        return(PI) 
        
      }
      
    }		
    
    entangle = function(PI){
      init_sum = sum(PI)
      for(i in 1:length(PI)){
        PI[i] = PI[i] - PI[i]/init_sum*E
      }
      return(PI)
    }    #unrealistic senario in which entanglment isn't effected by subpopulation size
    
    for(i in 1:(ncol(P))){
      
      P0 = P[,i]
      K = Ka[i]
      
      TP1 = sum(P0)
      
      Kef = Keffect(K,ev1,TP1)
      LK = L*Kef #apply effect of carrying capacity to Leslie matrix
      PI = LK%*%P0	 #population at t = 1 
      
      if(H > 0 && H_delay < 1){
        PI = apply_H(PI)
      } #apply hunting if delay less than 1 cycle is specified
      
      if(E>0 && sum(PI)<E){
        PI = entangle(PI)
      }
      
      PI = PI + (rbeta(PI, 5, 5)*2-1)*PI*R # +- random value within 5% of each age group based on a beta distribution a = b = 5
      PI = check(PI)
      
      P_out[,i] = PI
      
    }
    
    migration = function(P_out){
      
      P_mod = matrix(0, nrow = nrow(P),ncol = ncol(P))
      
      sum_P = c()
      for(i in 1:ncol(P_out)){
        sum_P[i] = sum(P_out[,i])
      } 
      
      sum_P_I = c()
      for(i in 1:length(sum_P)){   
        sum_P_I[i] = sum((M*sum_P)[,i])
      }
      
      for(i in 1:ncol(P_out)){
        current_str = P_out[,i]
        for( j in i:nrow(P_out)){
          P_out[j,i] = sum_P_I[i]*current_str[j]/sum_P[i]		
        } 
      }
      
      return(P_out)
      
    } #migration function
    
    if(sum(M)>0){
      P_out = migration(P_out)
    } #apply migration
    
    if(t==1){
      
      return(P_out) 
      
    }else{			
      
      for(i in 2:t){
        
        for(j in 1:(ncol(P))){
          
          PI = P_out[,j]
          K = Ka[j]
          
          RM = ((rbeta(PI, 5, 5)*2-1)*PI*R) # +- random value within 5% of each age group based on a beta distribution a = b = 5
          
          TP = sum(PI) 
          Kef = Keffect(K,ev1,TP) #calculate effect of carrying capacity
          LK = L*Kef #apply effect of carrying capacity to Leslie matrix
          PI = LK%*%PI
          
          if(H>0 && i>H_delay){
            PI = apply_H(PI)
          }	
          
          if(E>0 && sum(PI)<E){
            PI = entangle(PI)
          }
          
          PI = PI + RM #apply stocasticity
          PI = check(PI)
          
          P_out[,j] = PI
          
        }
        
        if(sum(M)>0){
          P_out = migration(P_out)
        } #apply migration
        
        if(P_D > 0 && runif(1, 0, 1)<P_D){
          impact = W_D*O_D
          P_out = P_out - P_out*impact
        } #determine if an outbreak of a disease occurs and apply, currently across the whole metapopulation
        
      }
      
      return(P_out)
      
    }
    
  }
  
  Pt = PGK(P)
  colnames(Pt) = colnames(P)
  return(Pt)
  
}

#################################