


library("tidyverse")
library("lubridate")
library("haven")
library("readxl")
library('willsutils') #devtools::install_github("wjchulme/willsutils")
import::from(magrittr, "%$%")



if("project_functions" %in% search()) detach(project_functions)

with((project_functions <- new.env(parent=as.environment("package:stats"))),
  {
    #list project-specific functions here:

    mad.p = function(x, na.rm=FALSE){
      mn <- mean(x, na.rm=na.rm)
      #len <- sum(!is.na(x))
      mean(abs(mn-x), na.rm=na.rm)
    }


    var.p = function(x, y = NULL, na.rm = FALSE){
      len<-sum(!is.na(x))
      var(x, y, na.rm = na.rm)*(len-1)/len
    }


    sd.p = function(x, na.rm = FALSE){
      len<-sum(!is.na(x))
      sd(x, na.rm = na.rm)*sqrt((len-1)/len)
    }


    qmultinom <- function(p, prob, is.prob.cumul = TRUE){

      if (is.null(names(prob))){
        warning('`prob` is not a named vector: assuming names are `seq_along(prob)`')
        names(prob) <- seq_along(prob)
      }

      if (length(unique(names(prob))) != length(prob) ){
        stop('`names(prob)` are not unique')
      }


      if (!is.prob.cumul)
        prob <- cumsum(prob/sum(prob))
      else if ( any(diff(prob)<0) | any(prob>1) | any(prob<0) | prob[length(prob)]!=1 )
        stop("if `is.prob.cumul == TRUE`, `prob` must be non-decreasing, 0-1 bounded, and the last value must be 1")

      as.numeric(unlist(lapply(p, function(p, prob){ names(prob)[min(which(p<=prob))]}, prob=prob)))
    }



    qmultinom_itpl <- function(p, prob, is.prob.cumul = TRUE){

      names(prob) <- seq_along(prob)

      if (!is.prob.cumul)
        prob <- cumsum(prob/sum(prob))
      else if ( any(diff(prob)<0) | any(prob>1) | any(prob<0) | prob[length(prob)]!=1 )
        stop("if `is.prob.cumul == TRUE`, `prob` must be non-decreasing, 0-1 bounded, and the last value must be 1")
      approx(x=c(0,prob), y=c(0,seq_along(prob)), xout=p, ties=min)$y
    }






    make.init.fixed.Q=function(n.states_cluster, out.rate){
      Quncorrected <- as.matrix(Matrix::bdiag(purrr::map(n.states_cluster, function(states){matrix(rep(out.rate, states^2), states)})))
      if(length(n.states_cluster)==1 & n.states_cluster[1]==1)
        Q0 <- matrix(1,nrow=1)
      else
        Q0 <- Quncorrected-diag(diag(Quncorrected))
      Q0-diag(rowSums(Q0))
    }


    hmodel2list <- function(hmodel, hmmdist = TRUE){

      if(!hmodel$hidden) stop("hmodel.object is not a Hidden Markov Model")

      .msm.LOOKUP <- data.frame(
        label = msm:::.msm.HMODELS,
        hmmname = c("hmmCat", "hmmIdent", "hmmUnif", "hmmNorm", "hmmLNorm", "hmmExp", "hmmGamma", "hmmWeibull", "hmmPois", "hmmBinom", "hmmBetaBinom",
                    "hmmTNorm", "hmmMETNorm", "hmmMEUnif", "hmmNBinom", "hmmBeta", "hmmT",
                    "hmmClmTNorm", "hmmClmTNorm7"),
        stringsAsFactors = FALSE
      )

      # makes a state-specific vector of parameters extracted from hmodel into a list of parameters (treating hmmCat as a special case)
      makeargslist <- function(params, label){
        # params = named vector of parameters for the distribution function
        #label = label (character) of the distribution function
        if(!(label %in% .msm.LOOKUP$label)) stop("Distribution ", label, " not currently supported for hmodel2list")
        if(label=="categorical")
          list(prob = params[names(params) %in% c("p", "p0", "pbase")], basecat = params[names(params)=="basecat"])
        else if(label=="identity")
          list(x = params[names(params) == "which"])
        else
          as.list(params)
      }

      labellist <- purrr::array_branch(hmodel$labels)
      paramlist <- split(hmodel$pars, list(hmodel$parout, hmodel$parstate))
      paramnestedlist <- mapply(makeargslist, paramlist, labellist, SIMPLIFY=FALSE, USE.NAMES=FALSE)
      distlist <- lapply(labellist, function(label){match.fun(.msm.LOOKUP$hmmname[.msm.LOOKUP$label==label])})

      if(hmodel$mv){
        hmmdistlist <- purrr::invoke_map(distlist, paramnestedlist)
        hmmdistnestedlist <- split(hmmdistlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
        msmlist <- lapply(hmmdistnestedlist, function(hmmdist){purrr::lift_dl(msm::hmmMV)(hmmdist)})

        if(hmmdist)
          msmlist
        else
          split(paramnestedlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
      } else {
        if(hmmdist)
          purrr::invoke_map(distlist, paramnestedlist)
        else
          split(paramnestedlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
      }

    }


    prob2thres <- function(prob, trunc){
      cumuprob <- cumsum(prob)
      msm::qtnorm(cumuprob[-length(cumuprob)], 0, 1, lower= -trunc, upper=trunc)
    }

    thres2prob <- function(thres, trunc){
      cumuprob<-c(0,msm::ptnorm(thres, 0, 1, lower= -trunc, upper=trunc), 1)
      diff(cumuprob)
    }

    thres2diff <- function(thres, trunc){
      diffs<-c(thres[1], diff(thres))
      names(diffs)<-c("thres1","diff12","diff23","diff34","diff45","diff56")
      diffs
    }

    prob2diff <- function(prob, trunc){
      cumuprob<-cumsum(prob)
      thres<-msm::qtnorm(cumuprob[-length(cumuprob)], 0, 1, lower= -trunc, upper=trunc)
      diffs<-c(thres[1], diff(thres))
      names(diffs)<-c("thres1","diff12","diff23","diff34","diff45","diff56")
      diffs
    }

    diff2prob <- function(diff, trunc){
      thres<-cumsum(diff)
      cumuprob<-c(0,msm::ptnorm(thres, 0, 1, lower= -trunc, upper=trunc), 1)
      diff(cumuprob)
    }



    cdf_distance2 <- function(x, y, norm){
      if(!is.numeric(x) | !is.numeric(y)) stop("x, y should be be numeric")

      cdfx <- ecdf(x)
      cdfy <- ecdf(y)

      if(!any(norm %in% c(1, 2, Inf))) stop("norm must be 1, 2 or Inf")
      domain <- sort(union(x, y))

      dif <- cdfx(domain)-cdfy(domain)
      fun <- approxfun(domain, dif, method='constant', yleft=0, yright=0, f=1)
      height <- fun(domain)[-length(domain)]

      if(norm==1L){
        height <- abs(height)
        width <- diff(domain)
        sum(height*width)
      } else
      if(norm==2L){
        height <- height^2
        width <- diff(domain)
        sqrt(sum(height*width))
      } else
      if(is.infinite(norm)){
        max(abs(height))
      }
    }



    wassdiscrete <- function(probsexp, probsobs, norm=1, support=seq_along(probsexp)){

      if(length(probsexp) != length(probsobs))stop("probsexp and probsobs should be same length")
      if(length(probsexp) != length(support))stop("support should be same length as probsexp and probsobs")
      #if(sum(probsexp) != 1 | sum(probsobs) != 1)stop("probsexp and probsobs should each sum to 1")
      if(!(norm>=1 | norm%%1 ==0 | is.infinite(norm))) stop("norm must be a positive integer or infinite")

      len <- length(support)
      cmlexp <- cumsum(probsexp)
      cmlobs <- cumsum(probsobs)

      difs <- abs(cmlexp[-len] - cmlobs[-len])

      if(is.infinite(norm)){
        D <- max(difs)
      } else {
        summand <- (difs*diff(support))^norm
        D <- sum(summand)^(1/norm)
      }

      if(any(is.na(probsexp)) | any(is.na(probsobs))){
        D <- NA
      } else if(!all.equal(sum(probsexp), 1) | !all.equal(sum(probsobs), 1)){
        D <- NA
      }

      D
    }



    stateprobMV <- function(Phi, hmodel, time, y, w){

    if(!is.matrix(y)) {
      y <- matrix(y, nrow=1)
    }
    if((length(time) != nrow(y))) stop("dimensions of time and y not compatible")
    if(w<1 | w>length(time)) stop("w must be between 1 and the number of observations available")

    J <- length(time)
    time <- time[seq_len(w)+(J-w)]
    y <- y[seq_len(w)+(J-w), ]


    if(w==1) {

      #This is a mixture model since we only have one time-point (one observation)

      S <- length(hmodel)
      weights <- rep(NA, S)
      for(s in seq_len(S)){
        weights[s] <- (hmodel[[s]]$Anxiety$pars[3:9][y[1]] +
                       hmodel[[s]]$Delusions$pars[3:9][y[2]] +
                       hmodel[[s]]$Depression$pars[3:9][y[3]] +
                       hmodel[[s]]$Hallucinations$pars[3:9][y[4]]
                       )/sum(!is.na(y))
      }
      p_tJ <- weights/sum(weights)
    } else {

      #we create a msm object for convenience so we can call the viterbi.msm
      #which, confusingly, also runs the forward-backward algorithm, which
      #calculates the probability of each state at each time point conditional on all the data
      #we take the probabilities for the last time-point (equivalent to the forward algorithm!)
      obj <- msm::msm(
        formula = y ~ time,
        #initprobs = initprobs,
        est.initprobs = TRUE,
        qmatrix = Phi,
        hmodel = hmodel,
        fixedpars = TRUE
      )

      p_t <- msm::viterbi.msm(obj)$pstate
      p_tJ <- p_t[nrow(p_t),]
    }
    names(p_tJ) <- rownames(Phi)
    p_tJ
  }






    # This takes a hmm (Phi, hmodel) and a history of symptoms (time, y),
    # and determines the vector of state probabilities at all observation times from w to max(t)
    # Is uses data from the previous w observations
    # output is a matrix of dimension w*k, where k is the number of states
    stateprobMV_FB <- function(Phi, hmodel, time, y, w){

      if(!is.matrix(y)) {
        y <- matrix(y, nrow=1)
      }
      if((length(time) != nrow(y))) stop("dimensions of time and y not compatible")
      if(w<1 | w>length(time)) stop("w must be between 1 and the number of observations available")

      J <- length(time)
      time <- time[seq_len(w)+(J-w)]
      y <- y[seq_len(w)+(J-w), ]


      if(w==1) {

        #This is a mixture model since we only have one time-point (one observation)

        S <- length(hmodel)
        weights <- rep(NA, S)
        for(s in seq_len(S)){
          weights[s] <- (hmodel[[s]]$Anxiety$pars[3:9][y[1]] +
                         hmodel[[s]]$Delusions$pars[3:9][y[2]] +
                         hmodel[[s]]$Depression$pars[3:9][y[3]] +
                         hmodel[[s]]$Hallucinations$pars[3:9][y[4]]
                         )/4
        }
        p_t <- t(weights/sum(weights))
        colnames(p_t) <-rownames(Phi)
        dplyr::bind_cols(tibble::tibble(time), tibble::as_tibble(p_t))
      } else {

        #we create a msm object for convenience so we can call the viterbi.msm
        #which, confusingly, also runs the forward-backward algorithm, which
        #calculates the probability of each state at each time point conditional on all the data
        #we take the probabilities for the last time-point (equivalent to the forward algorithm!)
        obj <- msm::msm(
          formula = y ~ time,
          #initprobs = initprobs,
          est.initprobs = TRUE,
          qmatrix = Phi,
          hmodel = hmodel,
          fixedpars = TRUE
        )

        vit <- msm::viterbi.msm(obj)
        time <- vit$time
        p_t <- vit$pstate
        colnames(p_t) <- rownames(Phi)
        dplyr::bind_cols(tibble::tibble(time), tibble::as_tibble(p_t))
      }

    }


    # This calculates the future state probabilities at u time units into the future,
    # given current state probabilities p and transition probabilities Phi
    # output is a vector of length k, where k is the number of states in the model
    futurestateprob <-
      Vectorize(
        function(p, Phi, u){
          p %*% expm::expm(u*Phi)
        },
        vectorize.args="u"
      )


    # This calculates the probabilities of ever passing through the states at u time units into the future,
    # given current state probabilities p and transition probabilities Phi
    # output is a vector of length k, where k is the number of states in the model
    futurepassprob <-
      Vectorize(
        function(p, Phi, u){
          msm::ppass.msm(qmatrix = Phi, start=p, tot=u)
        },
        vectorize.args="u"
      )


    # calculates the point in time u such that the probability of being in state k between now and u exceeds threshold tau
    # loops through
    select_next_time_incr <- function(p, Phi, tau, k, increment=30/24/60, max=100){
      u <- 0
      p_u <- rep(0, length(p))
      while( sum(p_u[k]) < tau | u < max) #must sample when combined risk of states k exceeds p=tau
      {
        u <- u+increment
        p_u <- futurepassprob(p, Phi, u)
      }

      u
    }

    # calculates the point in time u such that the probability of being in state k between now and u exceeds threshold tau
    # uses optimisation routine
    select_next_time <- function(p, Phi, tau, k, tol=0.000001){
      u <- tol

      optfunc <- function(u){
        (sum(futurepassprob(p, Phi, u)[k])-tau)^2
      }
      optim(u, optfunc, method="L-BFGS-B", lower=0, control=list(pgtol=tol))$par

    }


    # moved time forward to next acceptable sampling time, if sampling time is out-of-hours
    # assumes t is measured in decimal days
    adjust_next_time <- function(t, u, limits=c(3/24,3), interval=c(9/24,21/24)){

      mingap <- limits[1]
      maxgap <- limits[2]
      earliest <- interval[1]
      latest <- interval[2]

      u <- min(max(u,mingap),maxgap)
      #if(u<mingap) u <- mingap
      #if(u>maxgap) u <- maxgap
      next_time <- t+u
      if(next_time%%1 < earliest) next_time <- floor(next_time)+earliest
      if(next_time%%1 > latest) next_time <- ceiling(next_time)+earliest

      next_time
    }






  }
)
attach(project_functions);rm(project_functions)

search()




#assumes ordinal values

#qmultinom(p=c(0,0.2,0.5,0.8,0.9,1), prob=cumsum(prop.table(table(factor(c(1,2,3,4,4), levels=1:5)))), is.prob.cumul=TRUE)
#qmultinom(p=c(0,0.2,0.5,0.8,0.9,1), prob=c("1"=0.2, "2"=0.4, "3"=0.5, "4"=1), is.prob.cumul=TRUE)
#qmultinom(p=c(0, 0.1, 0.2, 0.4, 0.4999, 0.5, 0.50001, 0.6, 1), prob=c("1"=0.2, "2"=0.4, "3"=0.5, "4"=1), is.prob.cumul=TRUE)
#qmultinom(p=c(0, 0.1, 0.2, 0.4, 0.4999, 0.5, 0.50001, 0.6, 0.9), prob=c("1"=0.2, "2"=0.4, "3"=0.5, "4"=0.91), is.prob.cumul=TRUE)
#qmultinom(p=c(0, 0.1, 0.2, 0.4, 0.4999, 0.5, 0.50001, 0.6, 0.9), prob=c(0, 0.1, 0.3,1), is.prob.cumul=TRUE)

# prob2thres(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025), 5)
# thres2prob(prob2thres(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025), 5), 5)
# thres2diff(prob2thres(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025), 5), 5)
# prob2diff(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025), 5)
# cumsum(thres2diff(prob2thres(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025), 5), 5))
#
# sum(c(0.4,  0.3, 0.1, 0.1, 0.05, 0.025,0.025))

