#################################
#run model 'reps' number of times saving the mean curve and quasi extinction risk

proportion = TRUE
Hc = c(rep(0.01,5),rep(0.1,6)) #chose the hunting vector
Kc=c(2380,1370,1155,2204,2509,1363,1319,894,1645,905,1062) #chose your K values
reps = 10 #number of iterations

t_tot_vals = matrix(dat = NA, nrow = (1),ncol = (t+1))#holding vector
labs = c(NA)

ptm = proc.time()#timer begin

for(i in 1:reps){
  
  Pt = PopMod(t = t, A = A, 
              Kc = Kc, 
              P = P, B = B, S = S, R = R, 
              Hc = Hc, 
              H_delay = H_delay, H_sk = H_sk, proportion = proportion, 
              P_D = P_D, O_D = O_D, W_D = W_D)

  Pt = data.frame(Pt)
  Pt = transpose(Pt)
  t_tot_vals = rbind(t_tot_vals,Pt)
  labs = c(labs,names(data.frame(P)))
} #cycle through reps, popualate tot_vals with populaitons

t_tot_vals = cbind(labs,t_tot_vals)

Means = matrix(dat = NA, nrow = (t+1),ncol = (ncol(P)))
CIs =  matrix(dat = NA, nrow = (t+1),ncol = (ncol(P)))

Pl_Means = c()
Pl_CIs =  c()
Cols = c()
yrs = c()

CI = function(P,rep){
  m = rowSums(P)/rep
  i96 = 1.96*( sqrt(rowSums((P-m)^2)/rep) )
  return(i96)
} #95% confidence intervals

for(i in 1:ncol(P)){
  cl = names(data.frame(P))[i]
  dat = data.frame(t_tot_vals[which(t_tot_vals[1]==cl),])
  Means[,i] = c(colSums(dat[,-1])/reps)
  Pl_Means = c(Pl_Means,Means[,i])
  
  dat = transpose(dat[,-1])
  CIs[,i] = CI(dat,reps)
  Pl_CIs = c(Pl_CIs,Means[,i])
  
  Cols = c(Cols, rep(cl,(t+1)))
  yrs = c(yrs,seq(0:t))

}

Means = data.frame(Means)
CIs = data.frame(CIs)
names(Means) = names(data.frame(P))
names(CIs ) = names(data.frame(P))
Years = seq(0,t)  
Means = cbind(Years,Means)
CIs = cbind(Years,CIs)
proc.time() - ptm

#################################
#Plot average output

Pl = data.frame(cbind(yrs,Cols,Pl_Means,Pl_CIs))

names(Pl) = c("Year","Colony","Mean","CI")

ggplot(Pl, aes(x = yrs, y = Pl_Means, colour = Colony)) +
  geom_line() +
  theme_bw() + labs(x="Year", y=expression("Population Size"))
