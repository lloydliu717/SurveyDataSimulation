library(sem)            # needed for latent variable analysis

### ################################################################################################
### change local directory
# setwd("RS_sims")

### ################################################################################################
# set a random seed so the data generated produces consistent results
#   if you want to run different versions of this simulation, comment
#   out the following line by placing a "#" at the start of the line

# my.new.seed <- .Random.seed[271]   ## 271 = median of odd primes in length(.Random.seed) 
# sink(file="time_stamps.txt",append=TRUE)
# cat("Random Seed:\n")
# cat("············\n")
# print(my.new.seed)
# cat("\n")
# sink()
# print(my.new.seed)
# set.seed(my.new.seed)



### ################################################################################################
### focal variables to be adjusted for the different runs
# FIXED  number of manifest variables
#        number of reverse coded items (1-4)
#        top or bottom if odd
#         ·· 1 rev coded:  1-top, 1-bottom
#         ·· 2 rev coded:  top, both, bottom
#         ·· 3 rev coded:  2-top, 2-bottom
#         ·· 4 rev coded:  3-top, 2-2, 3-bottom
#            COUNT MULTIPIER = 10x
# FIXED  granularity
#        percent contamination
# FIXED  sample size
#        range for lambdas (∆=0.15 to set bottom-value of 0.40)
#        maximum lambda    (0.85,0.75,0.65)
#         ·· {0.4,0.85}, {0.5,0.85}, {0.6,0.85}, {0.7,0.85}
#         ·· {0.4,0.75}, {0.5,0.75}, {0.6,0.75}
#         ·· {0.4,0.65}, {0.5,0.65}
#            COUNT MULTIPIER = 9x
n.man <- 8
n.samp <- 400
gran <- 7

rev.item.list <- list(c(2),c(6),c(2,3),c(2,6),c(6,7),c(2,3,6),c(2,6,7),
                      c(2,3,4,6),c(2,3,6,7),c(2,6,7,8))


### ################################################################################################
### set the model
tmp.text <- paste("V",1:n.man,sep="")
tmp.text <- paste("   FAC1 --> ",tmp.text,sep="")
tmp.text <- paste(tmp.text,",   lambda",sep="")
tmp.text <- paste(tmp.text,1:n.man,sep="")
tmp.text <- paste(tmp.text,",  NA\n",sep="")
tmp.text <- paste(tmp.text,collapse="")
tmp.text <- paste("   FAC1 <-> FAC1, NA,         1\n",tmp.text,sep="",collapse="")
cat(tmp.text)
model.sem <- specifyModel(text=tmp.text,quiet=TRUE)

### ################################################################################################
### create the alternative models (different anchors vs. standardized latent variable)
tmp.text <- sub("NA,         1","xi,       NA",tmp.text)
cat(tmp.text)
for(i in 1:8) {
   tmp.text.alt <- sub(paste(c("lambda",i),collapse=""),"xyz",tmp.text)   
   tmp.text.alt <- sub("xyz,  NA","NA,  1",tmp.text.alt)   
   tmp.text.eval <- "model.sem.xyz <- specifyModel(text=tmp.text.alt,quiet=TRUE)"
   tmp.text.eval <- sub("xyz",i,tmp.text.eval)
   eval(parse(text=tmp.text.eval))
   }
rm(tmp.text,tmp.text.alt,tmp.text.eval,i)


### ################################################################################################
### initialize an empty data frame
var.names <- c("row.ct","sim.ct","n.man","gran","n.samp",
               "per.con","l.max","l.gap","n.rev","type.rev",
               "F.orig","chisq.orig","chisq.H0.orig",
               "CFI.orig","NNFI.orig","RMSEA.orig","SRMR.orig",
               paste(paste("st.lam.",1:n.man,sep=""),".orig",sep=""),
               "converge.type","F.min","chisq","df","chisq.H0","df.H0",
               "CFI","NNFI","RMSEA","SRMR",
               paste("st.lam.",1:n.man,sep=""))
n.rows <- 51*10*9*15

MDF <- data.frame(array(NA,dim=c(n.rows,length(var.names))))
names(MDF) <- var.names
rm(var.names)

### ################################################################################################
### initialize an empty matrix for the Lagrange multipliers
MI.mat <- array(NA,dim=c(n.man,n.man,n.rows))


########################################################################################
### establish the list of bad data runs
bad.list <- list(1)
bad.ctr <- 1

ctr <- 0
taus <- c(-1.5,-0.9,-0.3,0.3,0.9,1.5)

########################################################################################
### NOTE:  for the tryCatch() functions, the following models are tested:
###        1. test the latent var sd = 1 model
###        2. if that model is problematic, try it again, but after sampling
###           80% of the data to obtain a convering model (from which parameter estimates
###           are borrowed as starting values)
###        3. if this doesn't work after 20 tries, try running thru the various anchor item
###           models
for(jj in 0:50) {
  RS.per <- 2*jj/100
  for(kk in 1:3) {
    lam.max <- {95-kk*10}/100
    lam.mins <- seq(0.4,lam.max-0.15,0.10)
    for(kk.b in 1:length(lam.mins)) {
      lam.min <- lam.mins[kk.b]
      lambda <- c(rep.int(lam.max,4),rep.int(lam.min,4))
      for(ll in 1:length(rev.item.list)) {
        rev.items <- rev.item.list[[ll]]
        n.rev <- length(rev.items)
        rev.type <- length(which(rev.items < 5)) - length(which(rev.items > 4))
        for(ii in 1:15) {
          ### set the principle counter ··············································
          ctr <- ctr + 1
          ### generate the randomly generated data set ·······························
          xi <- rnorm(n.samp,mean=0,sd=1)
          delta <- array(rnorm(n.samp*n.man,mean=0,sd=1),dim=c(n.samp,n.man))
          x.mat <- xi %*% t(lambda) + delta %*% diag(sqrt(1-lambda^2))
          XX.mat <- array(0,dim=c(n.samp,n.man))
          for(i in 1:length(taus)) { XX.mat <- XX.mat + {taus[i] < x.mat}*1 }
          XX.mat <- as.data.frame(XX.mat)
          sem.out.orig <- sem(model.sem, cov(XX.mat), dim(XX.mat)[[1]],maxiter=10000)
          sem.stdcoef.orig <- stdCoef(sem.out.orig)
          sem.modin.orig <- modIndices(sem.out.orig)
          sem.summary.orig <- summary(sem.out.orig,fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
          XX.mat.r <- XX.mat
          tmp.val <- round(RS.per*n.samp)+1
          if(tmp.val < n.samp) {
            t.i <- tmp.val:n.samp   ## t.i = tmp.inds
            for(i in 1:n.rev) XX.mat.r[t.i,rev.items[i]] <- 6 - XX.mat.r[t.i,rev.items[i]]
            }
          ## set the convergence FLAG
          converge.val <- -1
          ## run/attempt-to-run the INITIAL model
          tryCatch({
              sem.output <- sem(model.sem, cov(XX.mat.r), dim(XX.mat.r)[[1]],maxiter=10000)
              do.chk <- {sem.output$criterion < 0}
              do.chk <- do.chk || sem.output$iter > 9999
              do.chk <- do.chk || length(which(diag(sem.output$vcov)<0)) > 0
              if(!do.chk && converge.val == -1) converge.val <- 0 },
            error = function(e) { do.chk <- TRUE; } )
          ## run/attempt-to-run the INITIAL model v/ selected starting values
          if(do.chk) {
            fifty <- 1
            tmp.chk <- TRUE
            while(tmp.chk && fifty <= 50) {
              XX.mat.80 <- XX.mat.r[sample(1:n.samp,0.8*n.samp),]
              tryCatch({
                  sem.out.tmp <- sem(model.sem, cov(XX.mat.80),dim(XX.mat.80)[[1]],maxiter=10000);
                  tmp.chk <- {sem.out.tmp$criterion < 0};
                  tmp.chk <- tmp.chk || sem.out.tmp$iter > 9999
                  tmp.chk <- tmp.chk || length(which(diag(sem.out.tmp$vcov)<0)) > 0},
                error = function(e) {tmp.chk <- TRUE;})
              fifty <- fifty + 1
              }
            if(!tmp.chk) {
              tmp.model <- model.sem
              tmp.model[2:17,3] <- summary(sem.out.tmp)$coeff[,1]
              sem.output <- sem(tmp.model, cov(XX.mat.r), dim(XX.mat.r)[[1]],maxiter=10000);
              do.chk <- {sem.output$criterion < 0};
              do.chk <- do.chk || sem.output$iter > 9999
              do.chk <- do.chk || length(which(diag(sem.output$vcov)<0)) > 0
              }
            }
          if(!do.chk && converge.val == -1) converge.val <- 800 + fifty - 1
          ## run/attempt-to-run the VARIATION models
          next.model <- 1
          while(do.chk && next.model <= n.man) {
            tmp.eval.txt <- "tmp.model.sem <- model.sem.xyz"
            eval(parse(text=sub("xyz",next.model,tmp.eval.txt)))
            tryCatch({
                sem.output <- sem(tmp.model.sem, cov(XX.mat.r),dim(XX.mat.r)[[1]],maxiter=10000);
                do.chk <- {sem.output$criterion < 0};
                do.chk <- do.chk || sem.output$iter > 9999
                do.chk <- do.chk || length(which(diag(sem.output$vcov)<0)) > 0},
              error = function(e) {do.chk <- TRUE;err.val<-"B";})
            next.model <- next.model + 1
            }
          if(!do.chk && converge.val == -1) converge.val <- next.model - 1
          ## populate the values into the Master Data Frame

          MDF$row.ct[ctr] <- ctr
          MDF$sim.ct[ctr] <- ii
          MDF$n.man[ctr] <- n.man
          MDF$gran[ctr] <- gran
          MDF$n.samp[ctr] <- n.samp
          MDF$per.con[ctr] <- RS.per
          MDF$l.max[ctr] <- lam.max
          MDF$l.gap[ctr] <- lam.max-lam.min
          MDF$n.rev[ctr] <- n.rev
          MDF$type.rev[ctr] <- rev.type

          MDF$F.orig[ctr] <- sem.out.orig$criterion
          MDF$chisq.orig[ctr] <- sem.summary.orig$chisq
          MDF$chisq.H0.orig[ctr] <- sem.summary.orig$chisqNull
          MDF$CFI.orig[ctr] <- sem.summary.orig$CFI
          MDF$NNFI.orig[ctr] <- sem.summary.orig$NNFI
          MDF$RMSEA.orig[ctr] <- sem.summary.orig$RMSEA[1]
          MDF$SRMR.orig[ctr] <- sem.summary.orig$SRMR
          names.1 <- paste(paste("st.lam.",1:n.man,sep=""),".orig",sep="")
          names.2 <- paste("lambda",1:n.man,"",sep="")
          MDF[ctr,names.1] <- sem.stdcoef.orig[names.2,2]

          if(!do.chk) {
            sem.summary <- summary(sem.output,fit.indices = c("RMSEA", "samp.nFI", "NNFI", "CFI", "SRMR"))
            sem.stdcoef <- stdCoef(sem.output)
            sem.modin <- modIndices(sem.output)
            MDF$converge.type[ctr] <- converge.val
            MDF$F.min[ctr] <- sem.output$criterion
            MDF$chisq[ctr] <- sem.summary$chisq
            MDF$df[ctr] <- sem.summary$df
            MDF$chisq.H0[ctr] <- sem.summary$chisqNull
            MDF$df.H0[ctr] <- sem.summary$dfNull
            MDF$CFI[ctr] <- sem.summary$CFI
            MDF$NNFI[ctr] <- sem.summary$NNFI
            MDF$RMSEA[ctr] <- sem.summary$RMSEA[1]
            MDF$SRMR[ctr] <- sem.summary$SRMR
            names.1 <- paste("st.lam.",1:n.man,sep="")
            names.2 <- paste("lambda",1:n.man,"",sep="")
            MDF[ctr,names.1] <- sem.stdcoef[names.2,2]

            MI.mat[,,ctr] <- sem.modin$mod.P[1:n.man,1:n.man]
            }
          if(do.chk || {1 < converge.val && converge.val < 800}) {
            # cat("\n• "); cat(RS.per); cat(" & "); cat(lam.min); cat(" •\n");
            bad.list[[bad.ctr]] <- c(RS.per,lam.max,lam.min,n.rev,rev.type)
            bad.list[[bad.ctr+1]] <- XX.mat.r
            bad.ctr <- bad.ctr + 2
            if(do.chk) MDF$converge.type[ctr] <- -9999
            }
          ### #################################################################################
          ### progress indicator
          if(ctr %% 20 == 0) {
           if(ctr %% 100 == 0) {
            if(ctr %% 500 == 0) {
             cat("| ");
             if(ctr %% 1000 == 0) {
              cat(ctr/1000);
              cat("k");
              }
             cat("\n");
             } else cat(":")
            } else cat("·")
           }
          ### #################################################################################
          ### close all for-loops
          }
        }
      }
    }
  }

rm(jj,kk,kk.b,ll,ii,i,ctr,next.model,converge.val,do.chk,tmp.chk)
rm(RS.per,tmp.val,t.i,lam.max,lam.mins,lam.min,lambda,rev.items,n.rev,rev.type)
rm(xi,delta,x.mat,XX.mat,XX.mat.r,fifty,XX.mat.80,names.1,names.2)
rm(sem.output,sem.out.tmp,tmp.model,tmp.eval.txt,sem.summary,sem.stdcoef,sem.modin)


########################################################################################
### update a specific model
###    code specific to following conditions
###    random seed = 123456789
###    jj = 0:50
###    ii = 1:15
# MDF[which(is.na(MDF$st.lam.4)),]

# XX.now <- bad.list[[2]]
# sem(model.sem.4, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);

# tmp.model <- model.sem
# tmp.model[c(2:4,6:17),3] <- summary(sem(model.sem.4, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))$coeff[2:16,1]
# tmp.out <- sem(tmp.model, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
# tmp.sum <- summary(tmp.out,fit.indices = c("RMSEA", "samp.nFI", "NNFI", "CFI", "SRMR"))
# tmp.stco <- stdCoef(tmp.out)
# tmp.mis <- modIndices(tmp.out)

# ctr <- which(is.na(MDF$st.lam.4))
# MDF$F.min[ctr] <- tmp.out$criterion
# MDF$chisq[ctr] <- tmp.sum$chisq
# MDF$df[ctr] <- tmp.sum$df
# MDF$chisq.H0[ctr] <- tmp.sum$chisqNull
# MDF$df.H0[ctr] <- tmp.sum$dfNull
# MDF$CFI[ctr] <- tmp.sum$CFI
# MDF$NNFI[ctr] <- tmp.sum$NNFI
# MDF$RMSEA[ctr] <- tmp.sum$RMSEA[1]
# MDF$SRMR[ctr] <- tmp.sum$SRMR
# names.1 <- paste("st.lam.",1:n.man,sep="")
# names.2 <- paste("lambda",1:n.man,"",sep="")
# MDF[ctr,names.1] <- tmp.stco[names.2,2]


########################################################################################
### fix the number of decimals in a particular variable
MDF$l.gap <- round(MDF$l.gap,2)
### fix the reversing type (pattern)
MDF$type.rev <- sign(MDF$type.rev)

########################################################################################
### calcuate the uncensored CFI
attach(MDF,warn.conflicts = FALSE)
uCFI <- {{chisq.H0-df.H0}-{chisq-df}}/{chisq.H0-df.H0}
uCFI.orig <- {{chisq.H0.orig-df.H0}-{chisq.orig-df}}/{chisq.H0.orig-df.H0}
detach(MDF)
MDF <- cbind(MDF,uCFI)
MDF <- cbind(MDF,uCFI.orig)
rm(uCFI,uCFI.orig)

########################################################################################
### calcuate Hoetler's N
###    N = 1 + chisq95/F
Hoetler95 <- 1 + qchisq(.95,MDF$df)/MDF$F.min
Hoetler95.orig <- 1 + qchisq(.95,MDF$df)/MDF$F.orig
MDF <- cbind(MDF,Hoetler95,Hoetler95.orig)
rm(Hoetler95,Hoetler95.orig)

# write.csv(MDF,file="tmp_MDF_161130.csv")


########################################################################################
### some explorations for how to process/extract-data-from the MI matrix
tmp.vec <- as.vector(lower.triangle(MI.mat[2:8,1:7,1342]))
tmp.vec <- sort(tmp.vec[which(tmp.vec > 0)],decreasing=TRUE)

t <- 1:28
t <- t^2/28
t <- sqrt(t)*sqrt(28)
plot(t,tmp.vec)
for(i in 9:15) {
   t <- i:28
   abline(lm(tmp.vec[t] ~ t))
   }

t <- 19:28
abline(lm(tmp.vec[t] ~ t))

del.vec <- tmp.vec[1:27]-tmp.vec[2:28]
del.vec <- sort(del.vec,decreasing=TRUE)
plot(del.vec)

t <- 18:27
lm.tmp <- lm(del.vec[t] ~ t)
plot(del.vec)
abline(lm.tmp)

z.vals <- del.vec[1:17] - predict(lm.tmp,data.frame(t=seq(1,17,1)))
z.vals <- {abs(z.vals)}/sd(lm.tmp$residuals)
min(which(z.vals < 16/3))

t <- 14:27
# t <- sqrt(t)*sqrt(27)
lm.tmp <- lm(del.vec[t] ~ t)

t <- 1:27
# t <- sqrt(t)*sqrt(27)
plot(t,del.vec)
abline(lm.tmp)

z.vals <- del.vec[1:13] - predict(lm.tmp,data.frame(t=seq(1,13,1)))
z.vals <- {abs(z.vals)}/sd(lm.tmp$residuals)
min(which(z.vals < 16/3))

t <- 5:27
# t <- sqrt(t)*sqrt(27)
lm.tmp <- lm(del.vec[t] ~ t)

z.vals <- del.vec[1:4] - predict(lm.tmp,data.frame(t=seq(1,4,1)))
z.vals <- {abs(z.vals)}/sd(lm.tmp$residuals)
min(which(z.vals < 16/3))

t <- 1:27
t.sq <- t^2/27
t.sr <- sqrt(t)*sqrt(27)
plot(t,del.vec,pch=20,col="blue")
points(t.sq,del.vec,pch=20,cex=0.6,col="darkred")
points(t.sr,del.vec,pch=20,cex=0.6,col="darkgreen")


########################################################################################
### some explorations for how to run models with different starting values
use.this.mod <- 2

XX.now <- bad.list[[use.this.mod]]
sem(model.sem, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.1, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.5, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.6, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.7, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.8, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);

sem(model.sem.2, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
stdCoef(sem(model.sem.2, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))
sem(model.sem.3, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
sem(model.sem.4, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);

tmp.model <- model.sem
tmp.model[c(2,4:17),3] <- summary(sem(model.sem.4, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))$coeff[2:16,1]
sem(tmp.model, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);
stdCoef(sem(tmp.model, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))

XX.now <- bad.list[[use.this.mod]][sample(1:400,0.8*400),]
sem(model.sem, cov(XX.now),dim(XX.now)[[1]],maxiter=10000);

tmp.model <- model.sem
tmp.model[2:17,3] <- summary(sem(model.sem, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))$coeff[,1]

XX.now <- bad.list[[use.this.mod]]
sem(tmp.model, cov(XX.now),dim(XX.now)[[1]],maxiter=10000)
stdCoef(sem(tmp.model, cov(XX.now),dim(XX.now)[[1]],maxiter=10000))




########################################################################################
### some initial model runs to assess impact
# pdf(file="n_rev_equal_1.pdf",width=4,height=4)
par(mfcol=c(1,2),mar={c(4, 4, 1, 1) + 0.1})
rev.list <- 1:-1
for(j in 1:3) {
# rev.list <- c(1,-1)
# for(j in 1:2) {
rev.num <- rev.list[j]
tmp.data <- MDF
# tmp.data <- na.omit(MDF)
per.sq <- tmp.data$per.con^2
per.cu <- tmp.data$per.con^3;    per.qu <- tmp.data$per.con^4;
per.qn <- tmp.data$per.con^5;    per.sx <- tmp.data$per.con^6;
tmp.data <- cbind(tmp.data,per.sq,per.cu,per.qu,per.qn,per.sx)

tmp.inds <- which(tmp.data$n.rev==1)
tmp.inds <- intersect(tmp.inds,which(tmp.data$type.rev==rev.num))
# tmp.inds <- intersect(tmp.inds,which(tmp.data$l.max==0.65))
# tmp.inds <- intersect(tmp.inds,which(tmp.data$l.gap==0.45))

tmp.data <- tmp.data[tmp.inds,]
# plot(tmp.data$per.con,tmp.data$uCFI,pch=20,cex=0.25)

# plot(jitter(tmp.data$per.con),tmp.data$uCFI,pch=20,cex=0.25)
color.list <- c("lightsalmon1","lightskyblue1","palegreen1","brown3")
color.list <- color.list[tmp.data$type.rev+2]
plot(jitter(tmp.data$per.con),tmp.data$uCFI,pch=20,cex=0.45,col=color.list,
     xlab="percent contamination",xlim=c(0,1),ylab="uCFI",bty="n",
     main=paste("reverse type =",rev.num))

tmp.valsA <- c(0.85,0.85,0.85,0.85,0.75,0.75,0.75,0.65,0.65)
tmp.valsB <- c(0.15,0.25,0.35,0.45,0.15,0.25,0.35,0.15,0.25)
tmp.cols <- c(paste("blue",4:1,sep=""),paste("red",4:2,sep=""),paste("green",4:3,sep=""))
x <- {0:100}/100
for(i in 1:9) {
   t.d <- tmp.data[which({tmp.data$l.max == tmp.valsA[i]}*{tmp.data$l.gap == tmp.valsB[i]} == 1),]
   lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu + per.qn + per.sx,data=t.d)
   coefs <- summary(lm.out)$coef
   y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + coefs[4,1]*x^3 + 
        coefs[5,1]*x^4 + coefs[6,1]*x^5 + coefs[7,1]*x^6
   lines(x,y,col=tmp.cols[i],lwd=1.25)
   }
}
# dev.off()



tmp.vals <- c(0.85,0.75,0.65)
x <- {0:100}/100
for(i in 1:3) {
   t.d <- tmp.data[which(tmp.data$l.max == tmp.vals[i]),]
   lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
   coefs <- summary(lm.out)$coef
   y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
        coefs[4,1]*x^3 + coefs[5,1]*x^4
   lines(x,y,col="blue3",lwd=2.5)
   }

tmp.vals <- c(0.15,0.25,0.35,0.45)
x <- {0:100}/100
for(i in 1:4) {
   t.d <- tmp.data[which(tmp.data$l.gap == tmp.vals[i]),]
   lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
   coefs <- summary(lm.out)$coef
   y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
        coefs[4,1]*x^3 + coefs[5,1]*x^4
   lines(x,y,col="red4",lwd=2.5)
   }


t.d <- tmp.data[which({tmp.data$type.rev == 1} * {tmp.data$l.max == 0.65}==1),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="green4",lwd=1.5)






t.d <- tmp.data[which(tmp.data$type.rev == 0),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="blue2")
t.d <- tmp.data[which({tmp.data$type.rev == 0} * {tmp.data$l.max == 0.85}==1),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="brown2")
t.d <- tmp.data[which({tmp.data$type.rev == 0} * {tmp.data$l.max == 0.75}==1),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="brown3")
t.d <- tmp.data[which({tmp.data$type.rev == 0} * {tmp.data$l.max == 0.65}==1),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="brown4")

t.d <- tmp.data[which(tmp.data$type.rev == 1),]
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="green3")

lm.out <- lm(uCFI ~ (per.con + per.sq + per.cu + per.qu + per.qn + per.sx) *
                    (l.max + l.gap),data=tmp.data)
lm.out <- lm(uCFI ~ (per.con + per.sq + per.cu + per.qu) *
                    (l.max + l.gap),data=tmp.data)
summary(lm.out)






color.list <- c("orange2","blue2","green3","brown3")
color.list <- color.list[10*tmp.data$l.gap-0.5]
# plot(jitter(tmp.data$per.con),tmp.data$uCFI.orig - tmp.data$uCFI,pch=20,cex=0.25,col=color.list)
plot(jitter(tmp.data$per.con),tmp.data$uCFI,pch=20,cex=0.25,col=color.list)







lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu + per.qn + per.sx +
                    l.max + l.gap + n.rev + as.factor(sign(tmp.data$type.rev)),data=tmp.data)
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu +
                    l.max + l.gap + n.rev + as.factor(sign(tmp.data$type.rev)),data=tmp.data)

tmp.inds <- which(tmp.data$n.rev==2)
tmp.data <- tmp.data[tmp.inds,]
plot(MDF$per.con[tmp.inds],MDF$uCFI.orig[tmp.inds] - MDF$uCFI[tmp.inds],pch=20,cex=0.4)
t.d <- tmp.data[which(tmp.data$l.max == 0.85),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y)
t.d <- tmp.data[which(tmp.data$l.max == 0.75),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y)
t.d <- tmp.data[which(tmp.data$l.max == 0.65),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.45),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red1")
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.35),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red2")
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.25),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red3")
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.15),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red4")

tmp.inds <- which(MDF$n.rev == 2)
tmp.data <- MDF[tmp.inds,]
tmp.inds <- which(tmp.data$type.rev == -2)
tmp.data <- tmp.data[tmp.inds,]
per.sq <- tmp.data$per.con^2
per.cu <- tmp.data$per.con^3
per.qu <- tmp.data$per.con^4
per.qn <- tmp.data$per.con^5
per.sx <- tmp.data$per.con^6
tmp.data <- cbind(tmp.data,per.sq,per.cu,per.qu,per.qn,per.sx)
color.list <- c("red","blue","green")
color.list <- color.list[10*tmp.data$l.max-5.5]
# plot(jitter(tmp.data$per.con),tmp.data$uCFI.orig - tmp.data$uCFI,pch=20,cex=0.25,col=color.list)
plot(jitter(tmp.data$per.con),tmp.data$uCFI,pch=20,cex=0.25,col=color.list)
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu + per.qn + per.sx +
                    l.max + l.gap,data=tmp.data)
summary(lm.out)
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu +
                    l.max + l.gap,data=tmp.data)
summary(lm.out)
lm.out <- lm(uCFI.orig - uCFI ~ (per.con + per.sq + per.cu + per.qu)*(l.max + l.gap) +
                    l.max * l.gap,data=tmp.data)
summary(lm.out)
lm.out <- lm(uCFI ~ (per.con + per.sq + per.cu + per.qu)*(l.max + l.gap) +
                    l.max * l.gap,data=tmp.data)
summary(lm.out)

t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.15) %in% which(tmp.data$l.max == 0.65),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="green1",lwd=1)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.25) %in% which(tmp.data$l.max == 0.65),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="green2",lwd=2)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.15) %in% which(tmp.data$l.max == 0.75),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="blue1",lwd=1)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.25) %in% which(tmp.data$l.max == 0.75),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="blue2",lwd=2)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.35) %in% which(tmp.data$l.max == 0.75),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="blue3",lwd=3)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.15) %in% which(tmp.data$l.max == 0.85),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red1",lwd=1)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.25) %in% which(tmp.data$l.max == 0.85),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red2",lwd=2)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.35) %in% which(tmp.data$l.max == 0.85),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red3",lwd=3)
t.d <- tmp.data[which(round(tmp.data$l.gap,2) == 0.45) %in% which(tmp.data$l.max == 0.85),]
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu ,data=t.d)
coefs <- summary(lm.out)$coef
x <- {0:100}/100
y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
     coefs[4,1]*x^3 + coefs[5,1]*x^4
lines(x,y,col="red4",lwd=4)



lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu +
                    l.max + l.gap + as.factor(sign(tmp.data$type.rev)),data=tmp.data[tmp.inds,])
lm.out <- lm(uCFI.orig - uCFI ~ per.con + per.sq + per.cu + per.qu +
                    l.max + l.gap,data=tmp.data[tmp.inds,])
summary(lm.out)


lm.out <- lm(NNFI ~ per.con + per.sq + per.cu + per.qu + per.qn + per.sx +
                    l.max + l.gap + n.rev + as.factor(sign(tmp.data$type.rev)),data=tmp.data)
lm.out <- lm(uCFI ~ per.con + per.sq + per.cu + per.qu +
                    l.max + l.gap + n.rev + as.factor(sign(tmp.data$type.rev)),data=tmp.data)
lm.out <- lm(NNFI ~ per.con + per.sq +
                    l.max + l.gap + n.rev + as.factor(sign(tmp.data$type.rev)),data=tmp.data)

# plot(MDF$per.con,abs(MDF$st.lam.2),pch=16,cex=0.4)
#
# plot(MDF[,1],MDF[,2],pch=16,cex=0.4)
# MDF <- data.frame(MDF)
# names(MDF) <- c("RS.per","NNFI")
# RS.sq <- MDF$RS.per^2
# RS.cu <- MDF$RS.per^3
# RS.qu <- MDF$RS.per^4
# RS.qn <- MDF$RS.per^5
# RS.sx <- MDF$RS.per^6
# lm.out <- lm(MDF$NNFI ~ MDF$RS.per + RS.sq + RS.cu + RS.qu + RS.qn + RS.sx)
# coefs <- summary(lm.out)$coef
# x <- {0:100}/100
# y <- coefs[1,1] + coefs[2,1]*x + coefs[3,1]*x^2 + 
#      coefs[4,1]*x^3 + coefs[5,1]*x^4 + coefs[6,1]*x^5 + coefs[7,1]*x^6
# lines(x,y)


      MI.mat[,,ctr] <- sem.modin$mod.P[1:pp,1:pp]






































########################################################################################
### calcuate the uncensored CFI
attach(data.mat,warn.conflicts = FALSE)
uCFI <- {{chisqNull-dfNull}-{chisq-df}}/{chisqNull-dfNull}
detach(data.mat)
data.mat <- cbind(data.mat,uCFI)

########################################################################################
### calcuate Hoetler's N
###    N = 1 + chisq95/F
Hoetler95 <- 1 + qchisq(.95,data.mat$df)/data.mat$F
data.mat <- cbind(data.mat,Hoetler95)

########################################################################################
### extract the Modification Indices from the 3D-matrix
### to a row of vectors for the lower-triangular values
LT <- lower.tri(MI.mat[,,1])
MIvec.mat <- array(NA,dim=c(num.rows,pp*{pp-1}/2))
for(ii in 1:num.rows) MIvec.mat[ii,] <- MI.mat[,,ii][LT]

file.append <- paste(c("run04_p",pp,"a"),collapse="")
write.files <- TRUE
if(write.files) {
   tmp.fn <- paste(c("data_mat_",file.append,".csv"),collapse="")
   write.csv(data.mat,file=tmp.fn)
   tmp.fn <- paste(c("lambda_mat_",file.append,".csv"),collapse="")
   write.csv(lambda.mat,file=tmp.fn)
   tmp.fn <- paste(c("modind_mat_",file.append,".csv"),collapse="")
   write.csv(MIvec.mat,file=tmp.fn)
   }

save.bad.files <- TRUE
if(save.bad.files && {bad.ctr > 2}) {
   bad.data.info <- NULL
   for(ii in 1:{{bad.ctr-1}/2}) {
      tmp.fn <- paste(c("BAD_DATA_SETS/run",file.append,"_",ifelse(ii<1000,"0",""),
                        ifelse(ii<100,"0",""),ifelse(ii<10,"0",""),ii,".csv"),collapse="")
      write.csv(bad.list[[2*ii]],file=tmp.fn)
      bad.data.info <- rbind(bad.data.info,bad.list[[2*ii-1]])
      }
   bad.data.info <- as.data.frame(bad.data.info)
   names(bad.data.info) <- names(data.mat)[1:7]
   tmp.fn <- paste(c("BAD_DATA_SETS/run",file.append,"_data_info.csv"),collapse="")
   write.csv(bad.data.info,file=tmp.fn)
   }




