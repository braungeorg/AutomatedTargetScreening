#Small and miscellaneous functions
###############################################################################
#' transformi
#' 
#' Used to derive logarithmic, normalized values for SAX analysis
#' 
#' @param y the array of values to be normalized
#' @return normalized array
#' @export
transformi = function(y){
  x=log10(y)
  (x-mean(x))/sd(x)
}

#'relative
#'
#'Transform data into % with highest value being 100
#'
#' @param y data to be transformed
#' @export
relative = function(y){
  Top = max(y,na.rm = T)
  x = y/Top*100
  return(x)
}

#' SAX_mindist
#' 
#' Uses a SAX distance matrix to calculate the MINDIST between two SAX strings
#' 
#' @param string1 First SAX string
#' @param string2 Second SAX string
#' @param Nr_of_data_points Number of data points which were summarized by the SAX strings
#' @export
SAX_mindist = function(string1,string2,Nr_of_data_points,SAX_reference_table){
  string1 = as.character(unlist(stringr::str_split(string1,"")))
  string2 = as.character(unlist(stringr::str_split(string2,"")))
  differences = rep(NA,length(string1))
  for(d in 1:length(string1)){
    differences[d] = as.numeric(SAX_reference_table[match(string1[d],letters),match(string2[d],letters)])^2
  }
  MINDIST = sqrt(Nr_of_data_points/length(string1)) * sqrt(sum(differences))
  return(MINDIST)
}

#'round_up
#'
#'Simply rounding up with one digit
#'
#' @param x number to be round up
#' @export
round_up = function(x){
  return(ceiling(x,digits=1))
}

#'letter to number
#'
#'Transform a letter into a number
#'
#' @param x the letter to be transformed
#' @export
letter2number <- function(x) {utf8ToInt(x) - utf8ToInt("a") + 1L}

#'SAX_consensus
#'
#'generate a consensus sequence from SAX strings
#'
#' @param vector_of_strings vector of SAX strings to build consensus
#' @export
SAX_consensus = function(vector_of_strings){
  if(length(vector_of_strings)==1){
    return(vector_of_strings)
  }
  split = stringr::str_split(vector_of_strings,"")
  comparison_table = as.data.frame(matrix(nrow=length(vector_of_strings),ncol=length(split[[1]])))
  for(i in 1:length(split)){
    comparison_table[i,] = unname(sapply(split[[i]],letter2number))
  }
  
  Consensus_numeric = rep(NA,length(split[[1]]))
  
  for(j in 1:length(comparison_table)){
    frequencies = as.data.frame(table(comparison_table[,j]))
    Consensus_numeric[j] = as.integer(sum(as.numeric(levels(frequencies$Var1))*as.numeric(frequencies$Freq))/sum(frequencies$Freq))
  }
  
  Consensus_letters = letters[Consensus_numeric]
  return(paste0(Consensus_letters,collapse = ""))
  
}

#'sum letters
#'
#'Counts the occurrence of the highest four letters in the selected alphabet
#'
#' @param x SAX sequence to analyze
#' @param alphabet_size the size of the alphabet selected for SAX
#' @export
sum_letters = function(x,alphabet_size){
  sums = rep(NA,length(x))
  for(i in 1:length(x)){
    sums[i] = sum(stringr::str_count(x[i],c(letters[(alphabet_size):(alphabet_size-3)])))
  }
  return(sums)
}

#'extract_first_two
#'
#'Get the two first values of a string
#'
#' @param x string to extract first two values
#' @export
extract_first_two = function(x){
  return(substr(x,1,2))
}

#'Remove_ion
#'
#'Remove additional information like _M+H separated by underscore
#'
#' @param x string like the compound name including ion information
#' @export
Remove_ion = function(x){
  split = unlist(str_split(x,"_"))
  new = paste0(split[1:2],collapse = "_")
  return(new)
}

#'errors
#'
#'function to handle errors from tryCatch functions within the code, especially for parallelization
#'
#'  @param e the error from the tryCatch function
#'  @export
errors = function(e){
  stopCluster(cl)
  results_path = paste0("./3_Results/",sample)
  if(!file.exists(results_path)){
    dir.create(results_path)
  }
  write(paste0(currentFunction,"\n",toString(e)), paste0(results_path,"/",Sys.Date(),"_ERROR_log.txt"))
  save.image(file=paste0(results_path,"/",Sys.Date(),"_ERROR.RData"))
  print(paste0(currentFunction,"\n",toString(e)))
  e <<- e
  stop(e)
}

#Major functions
################################################################################

#'calc_shift
#'
#'calculates the shift of expected and measured retention time of internal standards in the sample
#'
#' @param sample_path path to the sample / file batch
#' @param References vector of names of reference files
#' @param IS target list of internal standards
#' @param ppm_val allowed variability of target m/z in ppm
#' @param method_time the time of the MS measurement
#' @param FullscanMS1 respective scan filter in MS1
#' @param alphabet_size the alphabet size used for SAX functions
#' @param zigzag_trigger_threshold defines the maximum fraction of the length of the EIC which can consist of zigzag patterns before triggering a more inclusive background calculation
#' @param extended_baseline_factor a ratio applied to the baseline to define values which are considered as baseline
#' @param maximum_nr_of_a within SAX the letter a defines values with the lowest values and a valid peak shall only have a distinct number of these low values to avoid too flat peaks
#' @param minimum_nr_of_high_intensity_letters dependent on the selected SAX alphabet, the SAX sequence of a valid peak should have at least the respective number of letters defining high intensities
#' @param minimum_peak_area the value of a peak area which needs to be exceeded to be a valid peak
#' @param maximum_peak_width the maximum allowed peak width in minutes
#' @param i index of the standard from IS
#' @param q index of the reference from References
#' @export
calc_shift = function(sample_path,References,IS,ppm_val,method_time,FullscanMS1,alphabet_size,zigzag_trigger_threshold,extended_baseline_factor,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,minimum_peak_area,maximum_peak_width,i,q){
  
  if (T%in%grepl(".mzML", References)) {
    EIC = RaMS::grabMSdata(files = paste0(sample_path, "/", References[q]),
                           grab_what = "EIC",
                           mz = IS$`m/z`[i],
                           ppm = ppm_val,
                           rtrange = NULL,
                           prefilter = -1)
    
    EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
    
  } else {
    Chroma <- rawrr::readChromatogram(
      rawfile = paste0(sample_path, "/", References[q]),
      mass = IS$`m/z`[i],
      tol = ppm_val * 2,
      filter = FullscanMS1,
      type = "xic"
    )
    
    EIC <- data.frame("rt" = Chroma[[1]]$times, "int" = Chroma[[1]]$intensities)
  }
  
  EIC$int[EIC$int==Inf] = 0
  
  Filter = EIC$rt>=(IS$`retention time`[i]-2)&EIC$rt<=(IS$`retention time`[i]+2)
  Filter_background = EIC$rt>(IS$`retention time`[i]+3) | EIC$rt<(IS$`retention time`[i]-3)
  
  EIC_background = dplyr::filter(EIC,Filter_background==F)
  EIC = dplyr::filter(EIC,Filter==T)
  
  bad_chroma = F
  
  if(length(EIC$int)<3|length(EIC_background$int)<3){
    bad_chroma = T
  }
  
  if(bad_chroma==F){
    plot(EIC$int~EIC$rt,type="p",main=paste0(IS$identity[i]))
    timesplit = ceiling(length(EIC$rt)/50)
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(v in 1:50){
      EICs[[v]] = EIC[(Start:(Start+timesplit)),]
      EICs[[v]] = EICs[[v]][!is.na(EICs[[v]]$rt)&!is.na(EICs[[v]]$int),]
      timerange = max(EICs[[v]]$rt)-min(EICs[[v]]$rt)
      bandwidth = 0.1*timerange # 1 Min has a bandwidth of 0.05 for ksmooth (empirical value)
      EICs_smooth[[v]] = as.data.frame(ksmooth(EICs[[v]]$rt,EICs[[v]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[v]] = EICs_smooth[[v]][!is.na(EICs_smooth[[v]]$x)&!is.na(EICs_smooth[[v]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    lines(EIC_smooth$rt,EIC_smooth$int,col="red",type = "p")
    
    #quantile as benchmark to determine zigzag
    threshold = as.numeric(quantile(EIC_background$int,probs=c(0.1)))
    for(we in 1:nrow(EIC_background)){
      if(we == 1){
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = "b"
        } else {
          sequenci = "a"
        }
      } else {
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = append(sequenci,"b")
        } else {
          sequenci = append(sequenci,"a")
        }
      }
    }
    sequencis = paste0(sequenci,collapse = "")
    
    #Smooth EIC_background
    EICs_background = list()
    EICs_background_smooth = list()
    timesplit = ceiling(length(EIC_background$rt)/3)
    Start = 1
    for(v in 1:50){
      EICs_background[[v]] = EIC_background[(Start:(Start+timesplit)),]
      EICs_background[[v]] = EICs_background[[v]][!is.na(EICs_background[[v]]$rt)&!is.na(EICs_background[[v]]$int),]
      #plot(EICs_background[[v]]$int~EICs_background[[v]]$rt,type="p")
      timerange = max(EICs_background[[v]]$rt)-min(EICs_background[[v]]$rt)
      bandwidth = 0.1*timerange # 1 Min has a bandwidth of 0.05 for ksmooth (empirical value)
      EICs_background_smooth[[v]] = as.data.frame(ksmooth(EICs_background[[v]]$rt,EICs_background[[v]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_background_smooth[[v]] = EICs_background_smooth[[v]][!is.na(EICs_background_smooth[[v]]$x)&!is.na(EICs_background_smooth[[v]]$y),]
      #lines(EICs_background_smooth[[v]]$x,EICs_background_smooth[[v]]$y,type="p",col="red")
      Start=Start+timesplit
      if(Start>length(EIC_background$rt)){
        break
      }
    }
    
    EIC_background_smooth = dplyr::bind_rows(EICs_background_smooth)
    colnames(EIC_background_smooth) = c("rt","int")
    
    #Check for "aba" and replace with "ccc"
    zigzag_positions = as.data.frame(stringr::str_locate_all(sequencis,"aba"))
    
    if(nrow(zigzag_positions)>0){
      for(zi in 1:nrow(zigzag_positions)){
        Start = zigzag_positions$start[zi]
        End = zigzag_positions$end[zi]
        sequenci[Start:End] = "c"
      }
    }
    
    abline(v=EIC_background$rt[sequenci=="c"],col=adjustcolor(col="lightgrey",alpha.f = 0.2))
    
    if((length(sequenci[sequenci=="c"])>(zigzag_trigger_threshold*length(sequenci)))&(max(EIC_background$int,na.rm = T)!=min(EIC_background$int,na.rm=T))){
      Background = as.numeric(quantile(EIC_background$int[sequenci=="c"],probs=c(0.8)))
    } else if(nrow(EIC_smooth[EIC_smooth$int!=0,])>0){
      Background = as.numeric(min(EIC_background_smooth$int[EIC_background_smooth$int!=0],na.rm = T))
    } else {
      Background = 0
    }
    
    abline(h=Background,col="blue",lwd=3)
    
  } else {
    EIC_smooth = EIC
    EIC_background = data.frame("rt"=c(1,2,3),"int"=c(0.000001,0.000001,0.000001))
    Background = 0
  }
  
  if(length(EIC$rt)<=1){
    EIC = data.frame("rt"=c(1,2,3),"int"=c(0,0,0))
  }
  
  LOD = as.numeric(3*Background)
  LOQ = as.numeric(10*Background)
  
  if(is.na(LOD)|is.nan(LOD)){
    bad_chroma = T
  }
  
  if(bad_chroma==F){
    abline(h=LOD,col="green")
    abline(h=LOQ,col="orange")
  }
  
  guggl = seewave::SAX(EIC$int,alphabet_size = 3,PAA_number = 10)
  
  if(!(T%in%grepl("c",guggl))){
    bad_chroma = F
  }
  
  diffl = abs(EIC_smooth$rt-(min(EIC_smooth$rt,na.rm = T)+0.1))
  density = round(match(min(diffl,na.rm = T),diffl),digits = 0)
  
  #density is a measure of # of datapoints per 0.1 minutes
  
  SAX_result = rep("a",length(EIC$int))
  SAX_result_smooth = rep("a",length(EIC_smooth$int))
  
  SAX_result_small = SAX_result
  
  SAX_result_shift = SAX_result
  
  if(bad_chroma==F){
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>(Background*extended_baseline_factor)){
        SAX_result_shift[l] = "b"
      } else {
        SAX_result_shift[l] = "a"
      }
    }
    
    if(grepl(paste0(letters,collapse = "|"),IS$ID[i])){
      for(l in 1:length(EIC$int)){
        if(EIC$int[l]>LOD){
          SAX_result_small[l] = "b"
        } else {
          SAX_result_small[l] = "a"
        }
      }
    } else {
      for(l in 1:length(EIC$int)){
        if(EIC$int[l]>LOQ){
          SAX_result_small[l] = "b"
        } else {
          SAX_result_small[l] = "a"
        }
      }
    }
    
    if(length(EIC$int)<3|var(EIC$int)==0){
      bad_chroma = T
    }
    
    if(!(T%in%grepl("b",SAX_result_small))){
      bad_chroma = T
    }
  }
  
  if(bad_chroma == F){
    Peaklist = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"Bg_Start"=NA,"Bg_End"=NA)
    
    if(T%in%grepl("b",SAX_result_small)){
      Peak = match(max(EIC_smooth$int,na.rm = T),EIC_smooth$int)
      Peak_rt = EIC_smooth$rt[Peak]
      Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
      Cutoffs = c()
      if(!is.na(Peak)&Peak_int>LOD){
        remaining_Peaks = T
      } else {
        remaining_Peaks = F
      }
      while(remaining_Peaks == T){
        Peak_start = NA
        Peak_end = NA
        Cutoff_r = NA
        Cutoff_l = NA
        if(Peak>1){
          Begin_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak-1):1]
          Points_index = (Peak-1):1
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_start = Points_index[c-counter]
                break
              }
              Peak_start = Points_index[[c]]+1
              break
            }
            if(Check[c]>=(Begin_int*0.7)){
              SAX_result_smooth[Points_index[c]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_start = Points_index[c-7]
              if(c==7){
                Peak_start = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_start = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_start = 1
                } else {
                  Peak_start = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_start = 1
        }
        
        if(Peak < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak+1):length(SAX_result_smooth)]
          Points_index = (Peak+1):length(SAX_result_smooth)
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_end = Points_index[c-counter]
                break
              }
              Peak_end = Points_index[[c]]-1
              break
            }
            if(Check[c]>=(0.7*End_int)){
              SAX_result_smooth[Points_index[c]] = "b"
              SAX_result_smooth[Points_index[c-1]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_end = Points_index[c-7]
              if(c==7){
                Peak_end = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_end = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_end = length(SAX_result_smooth)
                } else {
                  Peak_end = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_end = length(SAX_result_smooth)
        }
        
        SAX_result_smooth[Peak] = "b"
        
        if(Peak_start > 1){
          Begin_int = EIC_smooth$int[Peak_start]
          Points_int = EIC_smooth$int[(Peak_start):1]
          Points_index = (Peak_start):1
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_l = Points_index[c-counter]
                break
              }
              Cutoff_l = Points_index[c-1]
              break
            }
            #print(paste("c: ",c,"Check: ",Check[c],"Check2: ",Check[c-1],"Decision1: ",Check[c]<(Check[c-1]),"counter: ",counter,"Decision2: ",Check[c]<threshold))
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_l = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              #Check for 60 % decrease over 0.5 minutes
              if(counter==(5*density+1)){
                Cutoff_l = Points_index[c-((5*density+1)-1)]
                break
              }
            }
          }
          if(is.na(Cutoff_l)){
            Cutoff_l = 1
          }
        } else {
          Cutoff_l = 1
        }
        
        if(Peak_end < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak_end]
          Points_int = EIC_smooth$int[Peak_end:length(SAX_result_smooth)]
          Points_index = Peak_end:length(SAX_result_smooth)
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_r = Points_index[c-counter]
                break
              }
              Cutoff_r = Points_index[c-1]
              break
            }
            #print(paste("c: ",c,"Check: ",Check[c],"Check2: ",Check[c-1],"Decision1: ",Check[c]<(Check[c-1]),"counter: ",counter,"Decision2: ",Check[c]<threshold))
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_r = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b" #-1 because Check = Points_index+1 in length
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              #Check for 60 % decrease over 0.5 minutes
              if(counter==(5*density+1)){
                Cutoff_r = Points_index[c-((5*density+1)-1)]
                break
              }
            }
          }
          if(is.na(Cutoff_r)){
            Cutoff_r = length(EIC_smooth$rt)
          }
        } else {
          Cutoff_r = length(EIC_smooth$rt)
        }
        
        Pairs = length(Cutoffs)/2
        if(Pairs>0){
          for(pa in 1:Pairs){
            pairis = data.frame("Start"=seq(1,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2),"End"=seq(2,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2))
            Range = Cutoffs[pairis$Start[pa]]:Cutoffs[pairis$End[pa]]
            if(Cutoff_l%in%Range){
              Cutoff_l = max(Range)
            }
            if(Cutoff_r%in%Range){
              Cutoff_r = min(Range)
            }
          }
        }
        
        Cutoffs = append(Cutoffs,c(Cutoff_l,Cutoff_r))
        
        Indizes = 1:length(EIC_smooth$int)
        
        Filter = SAX_result_smooth!="b"
        
        Indizes = Indizes[Filter]
        
        Peak = Indizes[match(max(EIC_smooth$int[SAX_result_smooth!="b"]),EIC_smooth$int[Filter])]
        
        Peak_rt = EIC_smooth$rt[Peak]
        
        Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
        
        if(!is.na(Peak)&Peak_int>LOD){
          remaining_Peaks = T
        } else {
          remaining_Peaks = F
        }
      }
      SAX_result_smooth[Cutoffs] = "a"
      #Not necessary because "a" is now considered as part of the peak
      #if(1%in%Cutoffs){
      #  SAX_result_smooth[1] = "b"
      #}
      #if(length(SAX_result_smooth)%in%Cutoffs){
      #  SAX_result_smooth[length(SAX_result_smooth)] = "b"
      #}
    }
    
    Start = 1
    
    while(Start<length(SAX_result_smooth)){
      peakbegin = match("b",SAX_result_smooth[Start:length(SAX_result_smooth)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result_smooth[peakbegin:length(SAX_result_smooth)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result_smooth)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      Starttime = EIC_smooth$rt[peakbegin]
      Endtime = EIC_smooth$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    SAX_result[EIC$int<=extended_baseline_factor*Background] = "a"
    
    Start = 1
    #If there are high intensities but no steady increase (like internal standards), widen the peak so that abbba not bbbbb
    while(Start<length(SAX_result)){
      peakbegin = match("b",SAX_result[Start:length(SAX_result)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result[peakbegin:length(SAX_result)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      if(peakbegin>1){
        if(EIC$int[peakbegin-1]<=(extended_baseline_factor*Background)){
          peakbegin = peakbegin - 1
        }
      }
      if(peakend<length(SAX_result)){
        if(EIC$int[peakend+1]<=(extended_baseline_factor*Background)){
          peakend = peakend + 1
        }
      }
      Starttime = EIC$rt[peakbegin]
      Endtime = EIC$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    #Since the widen-script is removing some cutoffs, make sure they are present in SAX_result
    Cutoffs_time = EIC_smooth$rt[Cutoffs]
    
    for(klo in 1:length(Cutoffs_time)){
      diff = abs(EIC$rt-Cutoffs_time[klo])
      time = EIC$rt[match(min(diff,na.rm = T),diff)]
      SAX_result[EIC$rt==time] = "a"
    }
    
    Start = 1
    End = length(SAX_result)
    Count = 1
    Already = 0
    
    while(5==5){
      not_a = match(letters[2:26],SAX_result[Start:End])+Already
      k = not_a[order(not_a,decreasing = F)[1]]
      k = k #dont use first a, since first and last b are already values which are not further decreasing
      #I disagree - Cutoffs are used as values which are relevant, but shared by multiple peaks, hence use first a
      if(!is.na(k)){
        k = k
      }
      l = match("a",SAX_result[k+1:End])+k-Start+Already+1
      if(!is.na(l)){
        l = l#l-1 #dont use last a, since first and last b are already values which are not further decreasing
        #I disagree - Cutoffs are used as values which are relevant, but shared by multiple peaks, hence use last a
      }
      #test for gaps
      m = k - 1
      if(is.na(m)){
        break
      }
      if(m==0){
        m = 1
      }
      if(abs(EIC$rt[k]-EIC$rt[m])>0.01){
        m = k
      }
      if(is.na(l)){
        if(mean(!is.na(not_a))!=0){ #removed Already==0
          l = End
          Peaklist$Start_RT[Count] = EIC$rt[m]
          Peaklist$End_RT[Count] = EIC$rt[l]
          level = 0.5
          #Peakl = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
          Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
          Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
          Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
          if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
          } else {
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
          }
          Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                              error=function(e){"too_small"})
          Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
          Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
          Peaklist$Height[Count] = max(EIC$int[m:l])
          #Peaklist$RT[Count] = EIC$rt[match(max(EIC$int[m:l]),EIC$int[m:l])+m-1]
          Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
          Peaklist$RT[Count] = median(Peaktop$rt)
          #Peaklist$RT[Count] = EIC_smooth$rt[match(max(EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]]),EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]])+match(T,EIC_smooth$rt>=EIC$rt[m])-1]
          id = m:l
          Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
          Count = Count+1
          Peaklist[Count,] = NA
          break
        } else{
          break
        }
      }
      Peaklist$Start_RT[Count] = EIC$rt[m]
      Peaklist$End_RT[Count] = EIC$rt[l]
      #Peakl = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
      level=0.5
      Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
      Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
      Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
      if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
      } else {
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
      }
      Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                          error=function(e){"too_small"})
      #Peaklist$RT[Count] = EIC$rt[match(max(EIC$int[m:l]),EIC$int[m:l])+m-1]
      Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
      Peaklist$RT[Count] = median(Peaktop$rt)
      #Peaklist$RT[Count] = EIC_smooth$rt[match(max(EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]]),EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]])+match(T,EIC_smooth$rt>=EIC$rt[m])-1]
      Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
      Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
      Peaklist$Height[Count] = max(EIC$int[m:l])
      id = m:l
      Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
      Start = l+1
      Count = Count+1
      Already = l
      Peaklist[Count,] = NA
    }
    Peaklist = Peaklist[!is.na(Peaklist$RT)&Peaklist$Sequence!="too_small",]
    
    if(nrow(Peaklist)==0){
      Peaklist[1,] = NA
    }
    
    Check_sequence = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters)&Peaklist$Area>minimum_peak_area
    
    if(grepl(paste0(letters,collapse = "|"),IS$ID[i])){
      Filter = Peaklist$Height >= LOD & !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence
    } else {
      Filter = Peaklist$Height >= LOQ & !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence
    }
    
    Peaklist_final = filter(Peaklist,Filter==T)
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"Bg_Start"=NA,"Bg_End"=NA)
    }
  } else {
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA)
  }
  Peaklist_final$Compound = IS$identity[i]
  Peaklist_final$mz = IS$`m/z`[i]
  if(!is.na(Peaklist_final$Height[1])){
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA)
      Peaklist_final$Compound = IS$identity[i]
      Peaklist_final$mz = IS$`m/z`[i]
    }
  }
  diff = abs(Peaklist_final$RT-IS$`retention time`[i])
  Peaklist_final = Peaklist_final[match(min(diff),diff),]
  
  if(bad_chroma==F&nrow(Peaklist_final[!is.na(Peaklist_final$Height),])>0){
    for(qlk in 1:nrow(Peaklist_final)){
      abline(v=Peaklist_final$Start_RT[qlk],col="purple")
      abline(v=Peaklist_final$End_RT[qlk],col="purple")
      Peak = EIC[EIC$rt>=Peaklist_final$Start_RT[qlk]&EIC$rt<=Peaklist_final$End_RT[qlk],]
      lines(Peak$rt,Peak$int,type="l",col="purple",lwd=2)
    }
  }
  
  return(Peaklist_final$RT-IS$`retention time`[i])
  
}

#'get_peaks_all
#'
#'Uses references to retrieve a list of possible peaks for the respective target
#'
#' @param sample_path path to the sample / file batch
#' @param sample the name of the sample / file batch
#' @param RT_range the retention time range defining the window checked for peaks and patterns in each direction
#' @param Multiple file defining peak patterns if multiple peaks occur in the m/z channel within RT_range
#' @param method_time the time of the MS measurement
#' @param FullscanMS1 the scan filter used for analysis of the main peaks
#' @param Scan_filters list of scan filters used for screening, mainly used for MS2 annotation
#' @param precursor_window the m/z window of precursor masses included when searching for MS2 scans in .mzML files
#' @param outer_search_window_multiple_peaks the maximum time in minutes before and after the aimed retention time, defining the search window if peak patterns are expected for a wider RT range
#' @param inner_search_window_multiple_peaks the narrow time in minutes before and after the aimed retention time, defining the search window if peak patterns are expected for a narrow RT range
#' @param outer_intensity_ratio_multiple_peaks the intensity ratio based on the maximum intensity within the outer RT range, used to filter valid peaks, usually a higher ratio
#' @param inner_intensity_ratio_multiple_peaks the intensity ratio based on the maximum intensity within the inner RT range, used to filter valid peaks, usually a lower ratio
#' @param References_refined the names of References, sorted by value if possible
#' @param is_centroided relevant for mzML files if scan mode is centroided ir not
#' @param compounds_shift_corrected the list of targets with retention times corrected
#' @param ppm_val allowed variability of target m/z in ppm
#' @param Times_shift_r the number of times the search window was extended to the right
#' @param Times_shift_l the number of times the search window was extended to the left
#' @param alphabet_size the alphabet size used for SAX functions
#' @param zigzag_trigger_threshold defines the maximum fraction of the length of the EIC which can consist of zigzag patterns before triggering a more inclusive background calculation
#' @param normal_background_quantile the quantile used to define background intensities from the background EIC
#' @param higher_background_quantile the quantile used to define background intensities from the background EIC if the chromatogram is showing quality flags 
#' @param extended_baseline_factor a ratio applied to the baseline to define values which are considered as baseline
#' @param maximum_nr_of_a within SAX the letter a defines values with the lowest values and a valid peak shall only have a distinct number of these low values to avoid too flat peaks
#' @param minimum_nr_of_high_intensity_letters dependent on the selected SAX alphabet, the SAX sequence of a valid peak should have at least the respective number of letters defining high intensities
#' @param minimum_peak_area the value of a peak area which needs to be exceeded to be a valid peak
#' @param maximum_peak_width the maximum allowed peak width in minutes
#' @param extended_baseline_factor a ratio applied to the baseline to define values which are considered as baseline
#' @param minimum_background_ratio the minimum ratio expected for a valid chromatogram of maximum intensity and the normal background quantile
#' @param width_smoothing a numeric value needed to define the density of smoothing for the EIC
#' @param maximum_nr_of_a within SAX the letter a defines values with the lowest values and a valid peak shall only have a distinct number of these low values to avoid too flat peaks
#' @param minimum_nr_of_high_intensity_letters dependent on the selected SAX alphabet, the SAX sequence of a valid peak should have at least the respective number of letters defining high intensities
#' @param minimum_peak_area the value of a peak area which needs to be exceeded to be a valid peak
#' @param maximum_peak_width the maximum allowed peak width in minutes
#' @param minimum_nr_of_datapoints_per_peak the minimal number of data points needed to define a valid peak
#' @param z index of the compound within compounds_shift_corrected to be analyzed
#' @param Reference index of the reference file to be used
#' @export
get_peaks_all = function(sample_path,sample,RT_range,Multiple,method_time,FullscanMS1,Scan_filters,precursor_window,outer_search_window_multiple_peaks,inner_search_window_multiple_peaks,outer_intensity_ratio_multiple_peaks,inner_intensity_ratio_multiple_peaks,References_refined,compounds_shift_corrected,ppm_val,Times_shift_r=0,Times_shift_l=0,use_shifting_search_window=F,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,extended_baseline_factor,minimum_background_ratio,width_smoothing,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,z,Reference){
  
  if(compounds_shift_corrected$MS1.MS2[z]=="MS2"){
    o = z-1
    while(grepl(paste0(letters,collapse="|"),compounds_shift_corrected$ID[o])){
      o = o-1
    }
    diff = abs(compounds_shift_corrected$`m/z`[o]-as.numeric(names(Scan_filters)))
    MS_filter = Scan_filters[[match(min(diff),diff)]]
  } else if(compounds_shift_corrected$MS1.MS2[z]=="MS1"&grepl(paste0(letters,collapse="|"),compounds_shift_corrected$ID[z])){
    o = z-1
    while(grepl(paste0(letters,collapse="|"),compounds_shift_corrected$ID[o])){
      o = o-1
    }
    if(compounds_shift_corrected$MS1.MS2[o]=="MS1"){
      MS_filter = FullscanMS1
    } else {
      diff = abs(compounds_shift_corrected$`m/z`[z]-as.numeric(names(Scan_filters)))
      MS_filter = Scan_filters[[match(min(diff),diff)]]
    }
  } else if(compounds_shift_corrected$MS1.MS2[z]=="MS1"){
    MS_filter = FullscanMS1
  } else {
    diff = abs(compounds_shift_corrected$`m/z`[z]-as.numeric(names(Scan_filters)))
    MS_filter = Scan_filters[[match(min(diff),diff)]]
  }
  
  if (T%in%grepl(".mzML", References_refined)) {
    if(compounds_shift_corrected$MS1.MS2[z]=="MS2"){
      EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",References_refined[Reference]),
                             grab_what = "MS2",
                             mz = compounds_shift_corrected$`m/z`[z],
                             ppm = ppm_val,
                             rtrange = NULL,
                             prefilter = -1)
      EIC = data.frame("rt"=EIC$MS2$rt,"int"=EIC$MS2$int,"mz"=EIC$MS2$fragmz,"prec"=EIC$MS2$premz)
      precursor = compounds_shift_corrected$`m/z`[o]
      EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
      EIC = EIC[(EIC$mz>=RaMS::pmppm(compounds_shift_corrected$`m/z`[z],ppm_val)[1]&EIC$mz<=RaMS::pmppm(compounds_shift_corrected$`m/z`[z],ppm_val)[2]),]
      EIC = EIC[,c(1:2)]
    } else {
      EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",References_refined[Reference]),
                             grab_what = "EIC",
                             mz = compounds_shift_corrected$`m/z`[z],
                             ppm = ppm_val,
                             rtrange = NULL,
                             prefilter = -1)
      EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
    }
  } else {
    Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",References_refined[Reference]),
                                     mass = compounds_shift_corrected$`m/z`[z],
                                     tol = ppm_val*2,
                                     filter = MS_filter,
                                     type = "xic")
    
    EIC <- data.frame("rt" = Chroma[[1]]$times, "int" = Chroma[[1]]$intensities)
  }
  
  Cutoff_shift_r = 0.15
  Cutoff_shift_l = -0.15
  
  Times_shift_r_old = Times_shift_r
  Times_shift_l_old = Times_shift_l
  
  EIC$int[EIC$int==Inf] = 0
  Filter_background = EIC$rt>(compounds_shift_corrected$`retention time`[z]+1+RT_range+Cutoff_shift_r*Times_shift_r) | EIC$rt<(compounds_shift_corrected$`retention time`[z]-1-RT_range+Cutoff_shift_l*Times_shift_l)
  Filter = EIC$rt < (compounds_shift_corrected$`retention time`[z]-RT_range+(Cutoff_shift_l*Times_shift_l)) | EIC$rt > (compounds_shift_corrected$`retention time`[z]+RT_range+(Cutoff_shift_r*Times_shift_r))
  ID = as.numeric(substr(compounds_shift_corrected$identity[z],1,4))
  if(is.na(ID)){
    ID = paste0("IS",substr(compounds_shift_corrected$identity[z],3,4))
  }
  if(ID%in%Multiple$ID){
    Filter = EIC$rt < (compounds_shift_corrected$`retention time`[z]-outer_search_window_multiple_peaks) | EIC$rt > (compounds_shift_corrected$`retention time`[z]+outer_search_window_multiple_peaks)
  }
  
  bad_chroma = F
  
  EIC_background = dplyr::filter(EIC,Filter_background==F)
  EIC = dplyr::filter(EIC,Filter==F)
  if(length(EIC$int)<3|length(EIC_background$int)<3){
    bad_chroma = T
  }
  
  if(bad_chroma==F){
    plot(EIC$int~EIC$rt,type="p",main=paste0(compounds_shift_corrected$identity[z]," r: ",Times_shift_r," l: ",Times_shift_l))
    timesplit = ceiling(length(EIC$rt)/width_smoothing)
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(i in 1:width_smoothing){
      EICs[[i]] = EIC[(Start:(Start+timesplit)),]
      EICs[[i]] = EICs[[i]][!is.na(EICs[[i]]$rt)&!is.na(EICs[[i]]$int),]
      timerange = max(EICs[[i]]$rt)-min(EICs[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_smooth[[i]] = as.data.frame(ksmooth(EICs[[i]]$rt,EICs[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[i]] = EICs_smooth[[i]][!is.na(EICs_smooth[[i]]$x)&!is.na(EICs_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    lines(EIC_smooth$rt,EIC_smooth$int,col="red",type = "p")
    
    threshold = as.numeric(quantile(EIC_background$int,probs=c(normal_background_quantile)))
    for(we in 1:nrow(EIC_background)){
      if(we == 1){
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = "b"
        } else {
          sequenci = "a"
        }
      } else {
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = append(sequenci,"b")
        } else {
          sequenci = append(sequenci,"a")
        }
      }
    }
    sequencis = paste0(sequenci,collapse = "")
    
    EICs_background = list()
    EICs_background_smooth = list()
    timesplit = ceiling(length(EIC_background$rt)/width_smoothing)
    Start = 1
    for(i in 1:width_smoothing){
      EICs_background[[i]] = EIC_background[(Start:(Start+timesplit)),]
      EICs_background[[i]] = EICs_background[[i]][!is.na(EICs_background[[i]]$rt)&!is.na(EICs_background[[i]]$int),]
      timerange = max(EICs_background[[i]]$rt)-min(EICs_background[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_background_smooth[[i]] = as.data.frame(ksmooth(EICs_background[[i]]$rt,EICs_background[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_background_smooth[[i]] = EICs_background_smooth[[i]][!is.na(EICs_background_smooth[[i]]$x)&!is.na(EICs_background_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC_background$rt)){
        break
      }
    }
    
    EIC_background_smooth = dplyr::bind_rows(EICs_background_smooth)
    colnames(EIC_background_smooth) = c("rt","int")
    
    zigzag_positions = as.data.frame(stringr::str_locate_all(sequencis,"aba"))
    
    if(nrow(zigzag_positions)>0){
      for(zi in 1:nrow(zigzag_positions)){
        Start = zigzag_positions$start[zi]
        End = zigzag_positions$end[zi]
        sequenci[Start:End] = "c"
      }
    }
    
    abline(v=EIC_background$rt[sequenci=="c"],col=adjustcolor(col="lightgrey",alpha.f = 0.2))
    
    if (T%in%grepl(".mzML", References_refined)){
      mean_Background = as.numeric(quantile(EIC_background$int,na.rm = T,probs = c(normal_background_quantile)))
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    } else {
      mean_Background = mean(EIC_background$int,na.rm = T)
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    }
    
    if((length(sequenci[sequenci=="c"])>(zigzag_trigger_threshold*length(sequenci)))&(max(EIC_background$int,na.rm = T)!=min(EIC_background$int,na.rm=T))){
      Q = as.numeric(quantile(EIC_background$int[sequenci=="c"],probs=c(higher_background_quantile)))
      Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
      Background = sd(Background_signals,na.rm=T)
      Baseline = mean(Background_signals)
    } else if(nrow(EIC_smooth[EIC_smooth$int!=0,])>0){
      if(Background_ratio<minimum_background_ratio){
        Q = as.numeric(quantile(EIC_background_smooth$int,probs=c(higher_background_quantile)))
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
          Background = sd(Background_signals,na.rm=T)
          Baseline = mean(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
        
      } else {
        Q = as.numeric(quantile(EIC_background_smooth$int,probs=c(normal_background_quantile)))
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<=Q]
          Baseline = mean(Background_signals,na.rm=T)
          Background = sd(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
      }
    } else {
      Background = 0
      Baseline = 0
    }
    
    abline(h=Baseline,col="blue",lwd=3)
    
  } else {
    EIC_smooth = EIC
    EIC_background = data.frame("rt"=c(1,2,3),"int"=c(0.000001,0.000001,0.000001))
    Background = 0
    Baseline = 0
  }
  
  if(length(EIC$rt)<=1){
    EIC = data.frame("rt"=c(1,2,3),"int"=c(0,0,0))
  }
  
  LOD = as.numeric(Baseline+3*Background)
  LOQ = as.numeric(Baseline+10*Background)
  
  if(is.na(LOD)|is.nan(LOD)){
    bad_chroma = T
  }
  
  if(bad_chroma==F){
    abline(h=LOD,col="green")
    abline(h=LOQ,col="orange")
  }
  
  guggl = seewave::SAX(EIC$int,alphabet_size = 3,PAA_number = 10)
  
  if(!(T%in%grepl("c",guggl))){
    bad_chroma = F
  }
  
  diffl = abs(EIC_smooth$rt-(min(EIC_smooth$rt,na.rm = T)+0.1))
  density = round(match(min(diffl,na.rm = T),diffl),digits = 0)
  
  SAX_result = rep("a",length(EIC$int))
  SAX_result_smooth = rep("a",length(EIC_smooth$int))
  
  SAX_result_small = SAX_result
  
  SAX_result_shift = SAX_result
  
  if(bad_chroma==F){
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>(Baseline*extended_baseline_factor)){
        SAX_result_shift[l] = "b"
      } else {
        SAX_result_shift[l] = "a"
      }
    }
    if(grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID[z])){
      for(l in 1:length(EIC$int)){
        if(EIC$int[l]>LOD){
          SAX_result_small[l] = "b"
        } else {
          SAX_result_small[l] = "a"
        }
      }
    } else {
      for(l in 1:length(EIC$int)){
        if(EIC$int[l]>LOQ){
          SAX_result_small[l] = "b"
        } else {
          SAX_result_small[l] = "a"
        }
      }
    }
    
    if(length(EIC$int)<3|var(EIC$int)==0){
      bad_chroma = T
    }
    
    if(!(T%in%grepl("b",SAX_result_small))){
      bad_chroma = T
    }
  }
  
  
  if(bad_chroma == F){
    Peaklist = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"level"=NA,"Concentration"=NA,"Bg_Start"=NA,"Bg_End"=NA)
    
    if(T%in%grepl("b",SAX_result_small)){
      Peak = match(max(EIC_smooth$int,na.rm = T),EIC_smooth$int)
      Peak_rt = EIC_smooth$rt[Peak]
      Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
      Cutoffs = c()
      if(!is.na(Peak)&Peak_int>LOD){
        remaining_Peaks = T
      } else {
        remaining_Peaks = F
      }
      while(remaining_Peaks == T){
        Peak_start = NA
        Peak_end = NA
        Cutoff_r = NA
        Cutoff_l = NA
        if(Peak>1){
          Begin_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak-1):1]
          Points_index = (Peak-1):1
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_start = Points_index[c-counter]
                break
              }
              Peak_start = Points_index[[c]]+1
              break
            }
            if(Check[c]>=(Begin_int*0.7)){
              SAX_result_smooth[Points_index[c]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_start = Points_index[c-7]
              if(c==7){
                Peak_start = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_start = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_start = 1
                } else {
                  Peak_start = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_start = 1
        }
        
        if(Peak < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak+1):length(SAX_result_smooth)]
          Points_index = (Peak+1):length(SAX_result_smooth)
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_end = Points_index[c-counter]
                break
              }
              Peak_end = Points_index[[c]]-1
              break
            }
            if(Check[c]>=(0.7*End_int)){
              SAX_result_smooth[Points_index[c]] = "b"
              SAX_result_smooth[Points_index[c-1]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_end = Points_index[c-7]
              if(c==7){
                Peak_end = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_end = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_end = length(SAX_result_smooth)
                } else {
                  Peak_end = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_end = length(SAX_result_smooth)
        }
        
        SAX_result_smooth[Peak] = "b"
        
        if(Peak_start > 1){
          Begin_int = EIC_smooth$int[Peak_start]
          Points_int = EIC_smooth$int[(Peak_start):1]
          Points_index = (Peak_start):1
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_l = Points_index[c-counter]
                break
              }
              Cutoff_l = Points_index[c-1]
              break
            }
            
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_l = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              
              if(counter==(5*density+1)){
                Cutoff_l = Points_index[c-((5*density+1)-1)]
                break
              }
            }
          }
          if(is.na(Cutoff_l)){
            Cutoff_l = 1
          }
        } else {
          Cutoff_l = 1
        }
        
        if(Peak_end < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak_end]
          Points_int = EIC_smooth$int[Peak_end:length(SAX_result_smooth)]
          Points_index = Peak_end:length(SAX_result_smooth)
          Check = Points_int
          counter = 0
          counter2 = 0
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_r = Points_index[c-counter]
                break
              }
              Cutoff_r = Points_index[c-1]
              break
            }
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_r = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              
              if(counter==(5*density+1)){
                Cutoff_r = Points_index[c-((5*density+1)-1)]
                break
              }
            }
          }
          if(is.na(Cutoff_r)){
            Cutoff_r = length(EIC_smooth$rt)
          }
        } else {
          Cutoff_r = length(EIC_smooth$rt)
        }
        
        Pairs = length(Cutoffs)/2
        if(Pairs>0){
          for(pa in 1:Pairs){
            pairis = data.frame("Start"=seq(1,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2),"End"=seq(2,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2))
            Range = Cutoffs[pairis$Start[pa]]:Cutoffs[pairis$End[pa]]
            if(Cutoff_l%in%Range){
              Cutoff_l = max(Range)
            }
            if(Cutoff_r%in%Range){
              Cutoff_r = min(Range)
            }
          }
        }
        
        Cutoffs = append(Cutoffs,c(Cutoff_l,Cutoff_r))
        
        Indizes = 1:length(EIC_smooth$int)
        
        Filter = SAX_result_smooth!="b"
        
        Indizes = Indizes[Filter]
        
        Peak = Indizes[match(max(EIC_smooth$int[SAX_result_smooth!="b"]),EIC_smooth$int[Filter])]
        
        Peak_rt = EIC_smooth$rt[Peak]
        
        Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
        
        if(!is.na(Peak)&Peak_int>LOD){
          remaining_Peaks = T
        } else {
          remaining_Peaks = F
        }
      }
      SAX_result_smooth[Cutoffs] = "a"
    }
    
    Start = 1
    
    while(Start<length(SAX_result_smooth)){
      peakbegin = match("b",SAX_result_smooth[Start:length(SAX_result_smooth)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result_smooth[peakbegin:length(SAX_result_smooth)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result_smooth)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      Starttime = EIC_smooth$rt[peakbegin]
      Endtime = EIC_smooth$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    SAX_result[EIC$int<=extended_baseline_factor*Baseline] = "a"
    
    Start = 1
    
    while(Start<length(SAX_result)){
      peakbegin = match("b",SAX_result[Start:length(SAX_result)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result[peakbegin:length(SAX_result)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      if(peakbegin>1){
        if(EIC$int[peakbegin-1]<=(extended_baseline_factor*Baseline)){
          peakbegin = peakbegin - 1
        }
      }
      if(peakend<length(SAX_result)){
        if(EIC$int[peakend+1]<=(extended_baseline_factor*Baseline)){
          peakend = peakend + 1
        }
      }
      Starttime = EIC$rt[peakbegin]
      Endtime = EIC$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    Cutoffs_time = EIC_smooth$rt[Cutoffs]
    
    for(klo in 1:length(Cutoffs_time)){
      diff = abs(EIC$rt-Cutoffs_time[klo])
      time = EIC$rt[match(min(diff,na.rm = T),diff)]
      SAX_result[EIC$rt==time] = "a"
    }
    
    Start = 1
    End = length(SAX_result)
    Count = 1
    Already = 0
    
    while(5==5){
      not_a = match(letters[2:26],SAX_result[Start:End])+Already
      k = not_a[order(not_a,decreasing = F)[1]]
      k = k 
      if(!is.na(k)){
        k = k
      }
      l = match("a",SAX_result[k+1:End])+k-Start+Already+1
      if(!is.na(l)){
        l = l
      }
      #test for gaps
      m = k - 1
      if(is.na(m)){
        break
      }
      if(m==0){
        m = 1
      }
      if(abs(EIC$rt[k]-EIC$rt[m])>0.01){
        m = k
      }

      if(is.na(l)){
        if(mean(!is.na(not_a))!=0){
          l = End
          level = 0.5
          Peaklist$Start_RT[Count] = EIC$rt[m]
          Peaklist$End_RT[Count] = EIC$rt[l]
          Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
          Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
          Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
          Peaklist$RT[Count] = median(Peaktop$rt)
          
          Peaklist$level[Count] = level
          
          Smoothie1 = Smoothie
          Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
          #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
          tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                                   resp = Smoothie1$int,
                                   bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
          
          RMSD_Ratio = Background/Baseline*100
          if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
            RMSD_Ratio = 0
          }
          if(RMSD_Ratio<10){
            Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
          } else {
            Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
          }
          
          if(is.na(Ratio)){
            Ratio = 0
          }
          
          #Constant = no real peak, too flat
          Constant_AIC = tcpl_fit$cnst_aic
          #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
          Gainloss_AIC = tcpl_fit$gnls_aic
          
          if(!is.na(Gainloss_AIC)){
            if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
              Gainloss_AIC = NA
            }
          }
          
          Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
          Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
          if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
          } else {
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
          }
          Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                              error=function(e){"too_small"})
          Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
          Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
          Peaklist$Height[Count] = max(EIC$int[m:l])
          
          id = m:l
          Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
          Count = Count+1
          Peaklist[Count,] = NA
          if(grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID[z])){
            if(Ratio<2|is.na(Gainloss_AIC)){
              Peaklist[(Count-1),] = NA
            }
          } else {
            if(Ratio<3|is.na(Gainloss_AIC)){
              Peaklist[(Count-1),] = NA
            }
          }
          break
        } else{
          break
        }
      }
      level = 0.5
      Peaklist$Start_RT[Count] = EIC$rt[m]
      Peaklist$End_RT[Count] = EIC$rt[l]
      Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
      Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
      Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
      Peaklist$RT[Count] = median(Peaktop$rt)
      Peaklist$level[Count] = level
      
      Smoothie1 = Smoothie
      Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
      #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
      tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                               resp = Smoothie1$int,
                               bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
      
      RMSD_Ratio = Background/Baseline*100
      if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
        RMSD_Ratio = 0
      }
      if(RMSD_Ratio<10){
        Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
      } else {
        Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
      }
      
      if(is.na(Ratio)){
        Ratio = 0
      }
      
      #Constant = no real peak, too flat
      Constant_AIC = tcpl_fit$cnst_aic
      #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
      Gainloss_AIC = tcpl_fit$gnls_aic
      
      if(!is.na(Gainloss_AIC)){
        if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
          Gainloss_AIC = NA
        }
      }
      
      Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
      Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
      if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
      } else {
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
      }
      Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                          error=function(e){"too_small"})
      Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
      Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
      Peaklist$Height[Count] = max(EIC$int[m:l])
      id = m:l
      Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
      Start = l+1
      Count = Count+1
      Already = l
      Peaklist[Count,] = NA
      if(grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID[z])){
        if(Ratio<2|is.na(Gainloss_AIC)){
          Peaklist[(Count-1),] = NA
        }
      } else {
        if(Ratio<3|is.na(Gainloss_AIC)){
          Peaklist[(Count-1),] = NA
        }
      }
    }
    Peaklist = Peaklist[!is.na(Peaklist$RT)&Peaklist$Sequence!="too_small",]
    
    if(nrow(Peaklist)==0){
      Peaklist[1,] = NA
    }
    
    Check_sequence_quan = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters)&Peaklist$Nr_of_Points>=(minimum_nr_of_datapoints_per_peak*2)&Peaklist$Area > minimum_peak_area
    Check_sequence_qual = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters)&Peaklist$Nr_of_Points>=minimum_nr_of_datapoints_per_peak&Peaklist$Area > (minimum_peak_area/2)
    
    if(grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID[z])){
      Filter = Peaklist$Height >= LOD & !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence_qual
    } else {
      Filter = Peaklist$Height >= LOQ & !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence_quan
    }
    
    Peaklist_final = filter(Peaklist,Filter==T)
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"level"=NA,"Concentration"=NA,"Bg_Start"=NA,"Bg_End"=NA)
    }
  } else {
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"level"=NA,"Concentration"=NA,"Bg_Start"=NA,"Bg_End"=NA)
  }
  Peaklist_final$Compound = compounds_shift_corrected$identity[z]
  Peaklist_final$mz = compounds_shift_corrected$`m/z`[z]
  if(!is.na(Peaklist_final$Height[1])){
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"level"=NA,"Concentration"=NA,"Bg_Start"=NA,"Bg_End"=NA)
      Peaklist_final$Compound = compounds_shift_corrected$identity[z]
      Peaklist_final$mz = compounds_shift_corrected$`m/z`[z]
    }
  }
  
  if(ID%in%Multiple$ID&nrow(Peaklist_final[!is.na(Peaklist_final$Height),])>1){
    #check if very small peak is included by mistake | threshold only applicable to Peaks with known patterns
    for(h in 1:nrow(Peaklist_final)){
      diff=abs(Peaklist_final$RT[h]-compounds_shift_corrected$`retention time`[z])
      if(diff>inner_search_window_multiple_peaks){
        if(Peaklist_final$Height[h]/max(Peaklist_final$Height,na.rm = T)<outer_intensity_ratio_multiple_peaks){
          Peaklist_final[h,] = NA
        } 
      } else {
        if(Peaklist_final$Height[h]/max(Peaklist_final$Height,na.rm = T)<inner_intensity_ratio_multiple_peaks){
          Peaklist_final[h,] = NA
        }
      }
    }
  }
  
  Peaklist_final = Peaklist_final[!is.na(Peaklist_final$Height),]
  if(nrow(Peaklist_final)==0){
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"level"=NA,"Concentration"=NA,"Bg_Start"=NA,"Bg_End"=NA)
    Peaklist_final$Compound = compounds_shift_corrected$identity[z]
    Peaklist_final$mz = compounds_shift_corrected$`m/z`[z]
  }
  
  if(bad_chroma==F&nrow(Peaklist_final[!is.na(Peaklist_final$Height),])>0){
    for(qlk in 1:nrow(Peaklist_final)){
      abline(v=Peaklist_final$Start_RT[qlk],col="purple")
      abline(v=Peaklist_final$End_RT[qlk],col="purple")
      Peak = EIC[EIC$rt>=Peaklist_final$Start_RT[qlk]&EIC$rt<=Peaklist_final$End_RT[qlk],]
      lines(Peak$rt,Peak$int,type="l",col="purple",lwd=2)
    }
  }
  
  if(!is.na(match("b",SAX_result))&!(ID%in%Multiple$ID)){
    if((SAX_result[2]=="b"&EIC$int[2]>(1.5*Baseline))|(EIC$int[(length(SAX_result)-1)]>(1.5*Baseline)&SAX_result[(length(SAX_result)-1)]=="b")){
      if(is.na(match("a",SAX_result))|(SAX_result[2]=="b"&EIC$int[2]>(1.5*Baseline))){
        Times_shift_l=Times_shift_l+1
      }
      
      if(is.na(match("a",SAX_result[length(SAX_result):1]))|(EIC$int[(length(SAX_result)-1)]>(1.5*Baseline)&SAX_result[(length(SAX_result)-1)]=="b")){
        Times_shift_r = Times_shift_r+1
      }
      
      if((Times_shift_r_old!=Times_shift_r|Times_shift_l_old!=Times_shift_l)&use_shifting_search_window==T){
        if(Times_shift_r<10&Times_shift_l<10){
          Peaklist_final = get_peaks_all(sample_path,sample,RT_range,Multiple,method_time,FullscanMS1,Scan_filters,precursor_window,outer_search_window_multiple_peaks,inner_search_window_multiple_peaks,outer_intensity_ratio_multiple_peaks,inner_intensity_ratio_multiple_peaks,References_refined,compounds_shift_corrected,ppm_val,Times_shift_r=Times_shift_r,Times_shift_l=Times_shift_l,use_shifting_search_window=T,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,extended_baseline_factor,minimum_background_ratio,width_smoothing,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,z,Reference)
        }
      }
    }
  }
  
  Peaklist_final$Bg_Start = min(EIC_background$rt,na.rm = T)
  Peaklist_final$Bg_End = max(EIC_background$rt,na.rm = T)
  
  if(grepl("CAL",References_refined[Reference])){
    Peaklist_final$density = density
    Concentration = gsub(".*CAL","CAL",References_refined[Reference])
    Concentration = gsub("p",".",Concentration)
    Concentration = as.numeric(substr(Concentration,5,8))
  } else {
    Peaklist_final$density = density
    Concentration = References_refined[Reference]
  }
  
  Peaklist_final$Concentration = Concentration
  
  return(Peaklist_final)
}

#'IS_analysis
#'
#'Checks for shift of internal standards per sample
#'
#' @param sample_path path to the sample / file batch
#' @param sample name of the sample / file batch
#' @param SAX_reference_table the SAX distance matrix used for calculating MINDIST
#' @param Solvent_Blank_ID the name ID of solvent blanks with no IS present
#' @param Multiple file defining peak patterns if multiple peaks occur in the m/z channel within RT_range
#' @param RT_range the retention time range defining the window checked for peaks and patterns in each direction
#' @param files names of the raw data files to be analyzed
#' @param Peaks results of peak annotation, data frame of peak characteristics from references
#' @param Peaks_list Peaks as list per reference
#' @param compounds_shift_corrected the list of targets with retention times corrected
#' @param FullscanMS1 the scan filter used for analysis of the main peaks
#' @param Scan_filters list of scan filters used for screening, mainly used for MS2 annotation
#' @param ppm_val allowed variability of target m/z in ppm
#' @param method_time the time of the MS measurement
#' @param is_centroided relevant for mzML files if scan is centroided or not
#' @param alphabet_size the alphabet size used for SAX functions
#' @param zigzag_trigger_threshold defines the maximum fraction of the length of the EIC which can consist of zigzag patterns before triggering a more inclusive background calculation
#' @param normal_background_quantile the quantile used to define background intensities from the background EIC
#' @param higher_background_quantile the quantile used to define background intensities from the background EIC if the chromatogram is showing quality flags 
#' @param minimum_background_ratio the minimum ratio expected for a valid chromatogram of maximum intensity and the normal background quantile
#' @param extended_baseline_factor a ratio applied to the baseline to define values which are considered as baseline
#' @param width_smoothing a numeric value needed to define the density of smoothing for the EIC
#' @param width_factor_background the factor applied to the peak width to define a background window
#' @param sample_search_window the search window in minutes added to both ends of the expected retention time used to screen for the peak of interest
#' @param maximum_nr_of_a within SAX the letter a defines values with the lowest values and a valid peak shall only have a distinct number of these low values to avoid too flat peaks
#' @param minimum_nr_of_high_intensity_letters dependent on the selected SAX alphabet, the SAX sequence of a valid peak should have at least the respective number of letters defining high intensities
#' @param minimum_peak_area the value of a peak area which needs to be exceeded to be a valid peak
#' @param maximum_peak_width the maximum allowed peak width in minutes
#' @param minimum_nr_of_datapoints_per_peak the minimal number of data points needed to define a valid peak
#' @param t index of the compound within Peaks
#' @param j index of the file to be analyzed in files
#' @export
IS_analysis = function(sample_path,sample,SAX_reference_table,Solvent_Blank_ID,Multiple,RT_range,files,Peaks,Peaks_list,compounds_shift_corrected,FullscanMS1,Scan_filters,ppm_val,method_time,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,minimum_background_ratio,extended_baseline_factor,width_smoothing,width_factor_background,sample_search_window,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,m,t,j){
  
  quan_peak = Peaks[t,]
  
  ID = as.numeric(substr(Peaks$Compound[t],1,4))
  if(is.na(ID)){ #Internal Standards
    ID = substr(Peaks$Compound[t],1,4)
  }
  
  if(is.null(quan_peak)){
    sample_peak = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA, "LOD" = NA, "LOQ" = NA)
    sample_peak$Compound = Peaks$Compound[t]
    sample_peak$mz = Peaks$mz[t]
    return(sample_peak)
    next
  }
  
  if(is.nan(quan_peak$RT)|is.na(quan_peak$RT)){
    sample_peak = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA, "LOD" = NA, "LOQ" = NA)
    sample_peak$Compound = Peaks$Compound[t]
    sample_peak$mz = Peaks$mz[t]
    return(sample_peak)
    next
  }
  
  if(compounds_shift_corrected$MS1.MS2[match(ID,compounds_shift_corrected$ID)]=="MS1"){
    MS_filter = FullscanMS1
  } else {
    diff = abs(Peaks$mz[t]-as.numeric(names(Scan_filters)))
    MS_filter = Scan_filters[[match(min(diff),diff)]]
  }
  
  if (T%in%grepl(".mzML", files)) {
    EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                           grab_what = "EIC",
                           mz = quan_peak$mz,
                           ppm = ppm_val,
                           rtrange = NULL,
                           prefilter = -1)
    EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
  } else {
    
    Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                     mass = quan_peak$mz,
                                     tol = ppm_val*2,
                                     filter = MS_filter,
                                     type = "xic")
    
    EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
  }
  
  EIC$int[EIC$int==Inf] = 0
  
  if(grepl(Solvent_Blank_ID,files[j])){
    maxis = quan_peak$RT
  } else {
    maxis = EIC$rt[which(EIC$int>0.5*max(EIC$int,na.rm = T))]
  }
  
  if(length(maxis)==0){
    maxis = quan_peak$RT
  }
  
  Background_window = quan_peak$Width*width_factor_background 
  
  if(Background_window<width_factor_background*sample_search_window){
    Background_window = width_factor_background*sample_search_window
  }
  
  maxi_frame = as.data.frame(matrix(nrow=nrow(EIC),ncol=length(maxis)))
  maxi_frame_background = as.data.frame(matrix(nrow=nrow(EIC),ncol=length(maxis)))
  
  for(m in 1:length(maxis)){
    maxi_frame[,m] = EIC$rt > (maxis[m]-quan_peak$Width-RT_range)&EIC$rt < (maxis[m]+quan_peak$Width+RT_range)
  }
  
  #Check that no peak is present at the expected time
  maxi_frame[,ncol(maxi_frame)+1] = EIC$rt > (quan_peak$RT-quan_peak$Width-RT_range)&EIC$rt < (quan_peak$RT+quan_peak$Width+RT_range)
  
  for(m in 1:length(maxis)){
    maxi_frame_background[,m] = EIC$rt > (maxis[m]-RT_range-Background_window) &EIC$rt < (maxis[m]+RT_range+Background_window)
  }
  
  maxi_frame_background[,ncol(maxi_frame_background)+1] = EIC$rt > (quan_peak$RT-RT_range-Background_window) &EIC$rt < (quan_peak$RT+RT_range+Background_window)
  
  Filter = rep(T,nrow(EIC))
  for(mi in 1:nrow(EIC)){
    if(T%in%maxi_frame[mi,]){
      Filter[mi] = F
    }
  }
  
  Filter_background = rep(F,nrow(EIC))
  
  for(mi in 1:nrow(EIC)){
    if(T%in%maxi_frame_background[mi,]){
      Filter_background[mi] = T
    }
  }
  
  bad_chroma = F
  
  EIC_background = dplyr::filter(EIC,Filter_background==T)
  EIC = dplyr::filter(EIC,Filter==F)
  if(length(EIC$rt)<3|length(EIC_background$rt)<3){
    bad_chroma = T
  }
  if(bad_chroma==F){
    plot(EIC$int~EIC$rt,type="p",main=quan_peak$Compound)
  }
  
  range = max(EIC$rt,na.rm = T)-min(EIC$rt,na.rm=T)
  
  timesplit = ceiling(length(EIC$rt)/width_smoothing)
  if(bad_chroma==F){
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(i in 1:width_smoothing){
      EICs[[i]] = EIC[(Start:(Start+timesplit)),]
      EICs[[i]] = EICs[[i]][!is.na(EICs[[i]]$rt)&!is.na(EICs[[i]]$int),]
      timerange = max(EICs[[i]]$rt)-min(EICs[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_smooth[[i]] = as.data.frame(ksmooth(EICs[[i]]$rt,EICs[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[i]] = EICs_smooth[[i]][!is.na(EICs_smooth[[i]]$x)&!is.na(EICs_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    lines(EIC_smooth$rt,EIC_smooth$int,col="red",type = "p")
    
    threshold = as.numeric(quantile(EIC_background$int,probs=c(normal_background_quantile)))
    for(we in 1:nrow(EIC_background)){
      if(we == 1){
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = "b"
        } else {
          sequenci = "a"
        }
      } else {
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = append(sequenci,"b")
        } else {
          sequenci = append(sequenci,"a")
        }
      }
    }
    sequencis = paste0(sequenci,collapse = "")
    
    EICs_background = list()
    EICs_background_smooth = list()
    timesplit = ceiling(length(EIC_background$rt)/width_smoothing)
    Start = 1
    for(i in 1:width_smoothing){
      EICs_background[[i]] = EIC_background[(Start:(Start+timesplit)),]
      EICs_background[[i]] = EICs_background[[i]][!is.na(EICs_background[[i]]$rt)&!is.na(EICs_background[[i]]$int),]
      timerange = max(EICs_background[[i]]$rt)-min(EICs_background[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_background_smooth[[i]] = as.data.frame(ksmooth(EICs_background[[i]]$rt,EICs_background[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_background_smooth[[i]] = EICs_background_smooth[[i]][!is.na(EICs_background_smooth[[i]]$x)&!is.na(EICs_background_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC_background$rt)){
        break
      }
    }
    
    EIC_background_smooth = dplyr::bind_rows(EICs_background_smooth)
    colnames(EIC_background_smooth) = c("rt","int")
    zigzag_positions = as.data.frame(stringr::str_locate_all(sequencis,"aba"))
    
    if(nrow(zigzag_positions)>0){
      for(zi in 1:nrow(zigzag_positions)){
        Start = zigzag_positions$start[zi]
        End = zigzag_positions$end[zi]
        sequenci[Start:End] = "c"
      }
    }
    
    if (T%in%grepl(".mzML", files)){
      mean_Background = as.numeric(quantile(EIC_background$int,na.rm = T,probs = c(normal_background_quantile)))
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    } else {
      mean_Background = mean(EIC_background$int,na.rm = T)
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    }
    
    if((length(sequenci[sequenci=="c"])>(zigzag_trigger_threshold*length(sequenci)))&(max(EIC_background$int,na.rm = T)!=min(EIC_background$int,na.rm=T))){
      Q = as.numeric(quantile(EIC_background$int[sequenci=="c"],probs=c(0.8)))
      Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
      Background = sd(Background_signals,na.rm=T)
      Baseline = mean(Background_signals)
    } else if(nrow(EIC_smooth[EIC_smooth$int!=0,])>0){
      if(Background_ratio<minimum_background_ratio){
        Q = as.numeric(quantile(EIC_background_smooth$int,probs=c(higher_background_quantile)))
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
          Background = sd(Background_signals,na.rm=T)
          Baseline = mean(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
        
      } else {
        Q = as.numeric(quantile(EIC_background_smooth$int,probs=c(normal_background_quantile)))
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<=Q]
          Baseline = mean(Background_signals,na.rm=T)
          Background = sd(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
      }
    } else {
      Background = 0
      Baseline = 0
    }
    
    nrow(EIC_smooth[EIC_smooth$int<=extended_baseline_factor*min(EIC_smooth$int),])>(zigzag_trigger_threshold*length(sequenci))
    
    abline(h=Background,col="blue",lwd=3)
    
  } else {
    EIC_smooth = EIC
    EIC_background = data.frame("rt"=c(1,2,3),"int"=c(0.000001,0.000001,0.000001))
    Background = 0
    Baseline = 0
  }
  
  if(length(EIC$rt)<=1){
    EIC = data.frame("rt"=c(1,2,3),"int"=c(0,0,0))
    Background = 0
    Baseline = 0
  }
  
  LOD = as.numeric(Baseline+3*Background)
  LOQ = as.numeric(Baseline+10*Background)
  if(bad_chroma==F){
    abline(h=LOD,col="green")
    abline(h=LOQ,col="orange")
  }
  
  if(is.na(LOD)|is.nan(LOD)){
    bad_chroma = T
  }
  
  guggl = seewave::SAX(EIC$int,alphabet_size = 3,PAA_number = 10)
  
  if(!(T%in%grepl("c",guggl))){
    bad_chroma = F
  }
  
  diffl = abs(EIC_smooth$rt-(min(EIC_smooth$rt,na.rm = T)+0.1))
  density = round(match(min(diffl,na.rm = T),diffl),digits = 0)
  
  SAX_result = rep("a",length(EIC$int))
  
  SAX_result_smooth = rep("a",length(EIC_smooth$int))
  
  SAX_result_small = SAX_result
  
  SAX_result_shift = SAX_result
  
  if(bad_chroma==F){
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>LOD){
        SAX_result_small[l] = "b"
      } else {
        SAX_result_small[l] = "a"
      }
    }
    
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>(Background*extended_baseline_factor)){
        SAX_result_shift[l] = "b"
      } else {
        SAX_result_shift[l] = "a"
      }
    }
    
    if(length(EIC$int)<3|var(EIC$int)==0){
      bad_chroma = T
    }
    
    if(!(T%in%grepl("b",SAX_result_small))){
      bad_chroma = T
    }
  }
  
  if(bad_chroma == F){
    Peaklist = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA)
    
    if(T%in%grepl("b",SAX_result_small)){
      Peak = match(max(EIC_smooth$int,na.rm = T),EIC_smooth$int)
      Peak_rt = EIC_smooth$rt[Peak]
      Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
      Cutoffs = c()
      if(!is.na(Peak)&Peak_int>LOQ){
        remaining_Peaks = T
      } else {
        remaining_Peaks = F
      }
      while(remaining_Peaks == T){
        Peak_start = NA
        Peak_end = NA
        Cutoff_r = NA
        Cutoff_l = NA
        if(Peak>1){
          Begin_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak-1):1]
          Points_index = (Peak-1):1
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_start = Points_index[c-counter]
                break
              }
              Peak_start = Points_index[[c]]+1
              break
            }
            if(Check[c]>=(Begin_int*0.7)){
              SAX_result_smooth[Points_index[c]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_start = Points_index[c-7]
              if(c==7){
                Peak_start = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_start = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_start = 1
                } else {
                  Peak_start = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_start = 1
        }
        
        if(Peak < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak+1):length(SAX_result_smooth)]
          Points_index = (Peak+1):length(SAX_result_smooth)
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_end = Points_index[c-counter]
                break
              }
              Peak_end = Points_index[[c]]-1
              break
            }
            if(Check[c]>=(0.7*End_int)){
              SAX_result_smooth[Points_index[c]] = "b"
              SAX_result_smooth[Points_index[c-1]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_end = Points_index[c-7]
              if(c==7){
                Peak_end = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_end = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_end = length(SAX_result_smooth)
                } else {
                  Peak_end = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_end = length(SAX_result_smooth)
        }
        
        SAX_result_smooth[Peak] = "b"
        
        if(Peak_start > 1){
          Begin_int = EIC_smooth$int[Peak_start]
          Points_int = EIC_smooth$int[(Peak_start):1]
          Points_index = (Peak_start):1
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_l = Points_index[c-counter]
                break
              }
              Cutoff_l = Points_index[c-1]
              break
            }
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_l = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              if(counter==(5*density+1)){
                Cutoff_l = Points_index[(c-((5*density+1)-1))]
                break
              }
            }
          }
          if(is.na(Cutoff_l)){
            Cutoff_l = 1
          }
        } else {
          Cutoff_l = 1
        }
        
        if(Peak_end < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak_end]
          Points_int = EIC_smooth$int[Peak_end:length(SAX_result_smooth)]
          Points_index = Peak_end:length(SAX_result_smooth)
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_r = Points_index[c-counter]
                break
              }
              Cutoff_r = Points_index[c-1]
              break
            }
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_r = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              if(counter==(5*density+1)){
                Cutoff_r = Points_index[(c-((5*density+1)-1))]
                break
              }
            }
          }
          if(is.na(Cutoff_r)){
            Cutoff_r = length(EIC_smooth$rt)
          }
        } else {
          Cutoff_r = length(EIC_smooth$rt)
        }
        
        Pairs = length(Cutoffs)/2
        if(Pairs>0){
          for(pa in 1:Pairs){
            pairis = data.frame("Start"=seq(1,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2),"End"=seq(2,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2))
            Range = Cutoffs[pairis$Start[pa]]:Cutoffs[pairis$End[pa]]
            if(Cutoff_l%in%Range){
              Cutoff_l = max(Range)
            }
            if(Cutoff_r%in%Range){
              Cutoff_r = min(Range)
            }
          }
        }
        
        Cutoffs = append(Cutoffs,c(Cutoff_l,Cutoff_r))
        
        Indizes = 1:length(EIC_smooth$int)
        
        Filter = SAX_result_smooth!="b"
        
        Indizes = Indizes[Filter]
        
        Peak = Indizes[match(max(EIC_smooth$int[SAX_result_smooth!="b"]),EIC_smooth$int[Filter])]
        
        Peak_rt = EIC_smooth$rt[Peak]
        
        Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
        
        if(!is.na(Peak)&Peak_int>LOQ){
          remaining_Peaks = T
        } else {
          remaining_Peaks = F
        }
      }
      SAX_result_smooth[Cutoffs] = "a"
      if(1%in%Cutoffs){
        SAX_result_smooth[1] = "b"
      }
      if(length(SAX_result_smooth)%in%Cutoffs){
        SAX_result_smooth[length(SAX_result_smooth)] = "b"
      }
    }
    
    Start = 1
    
    while(Start<length(SAX_result_smooth)){
      peakbegin = match("b",SAX_result_smooth[Start:length(SAX_result_smooth)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result_smooth[peakbegin:length(SAX_result_smooth)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result_smooth)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      Starttime = EIC_smooth$rt[peakbegin]
      Endtime = EIC_smooth$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    SAX_result[EIC$int<=extended_baseline_factor*Baseline] = "a"
    
    Start = 1
    while(Start<length(SAX_result)){
      peakbegin = match("b",SAX_result[Start:length(SAX_result)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result[peakbegin:length(SAX_result)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      if(peakbegin>1){
        if(EIC$int[peakbegin-1]<=(extended_baseline_factor*Baseline)){
          peakbegin = peakbegin - 1
        }
      }
      if(peakend<length(SAX_result)){
        if(EIC$int[peakend+1]<=(extended_baseline_factor*Baseline)){
          peakend = peakend + 1
        }
      }
      Starttime = EIC$rt[peakbegin]
      Endtime = EIC$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    Cutoffs_time = EIC_smooth$rt[Cutoffs]
    
    for(klo in 1:length(Cutoffs_time)){
      diff = abs(EIC$rt-Cutoffs_time[klo])
      time = EIC$rt[match(min(diff,na.rm = T),diff)]
      SAX_result[EIC$rt==time] = "a"
    }
    
    Start = 1
    End = length(SAX_result)
    Count = 1
    Already = 0
    
    while(5==5){
      not_a = match(letters[2:26],SAX_result[Start:End])+Already
      k = not_a[order(not_a,decreasing = F)[1]]
      k = k
      if(!is.na(k)){
        k = k
      }
      l = match("a",SAX_result[k+1:End])+k-Start+Already+1
      if(!is.na(l)){
        l = l
      }
      #test for gaps
      m = k - 1
      if(is.na(m)){
        break
      }
      if(m==0){
        m = 1
      }
      if(abs(EIC$rt[k]-EIC$rt[m])>0.01){
        m = k
      }
      if(is.na(l)){
        if(mean(!is.na(not_a))!=0){
          level = quan_peak$level
          l = End
          Peaklist$Start_RT[Count] = EIC$rt[m]
          Peaklist$End_RT[Count] = EIC$rt[l]
          Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
          Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
          Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
          Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
          if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
          } else {
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
          }
          Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                              error=function(e){"too_small"})
          Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
          Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
          Peaklist$Height[Count] = max(EIC$int[m:l])
          Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
          Peaklist$RT[Count] = median(Peaktop$rt)
          
          Smoothie1 = Smoothie
          Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
          #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
          tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                                   resp = Smoothie1$int,
                                   bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
          
          RMSD_Ratio = Background/Baseline*100
          if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
            RMSD_Ratio = 0
          }
          if(RMSD_Ratio<10){
            Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
          } else {
            Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
          }
          
          if(is.na(Ratio)){
            Ratio = 0
          }
          Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
          
          if(grepl("GC",sample)){
            if(Differenci<=0.01&Ratio<3){
              Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
            }
          } else {
            if(Differenci<=0.05&Ratio<3){
              Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
            }
          }
          
          #Constant = no real peak, too flat
          Constant_AIC = tcpl_fit$cnst_aic
          #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
          Gainloss_AIC = tcpl_fit$gnls_aic
          
          if(!is.na(Gainloss_AIC)){
            if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
              Gainloss_AIC = NA
            }
          }
          id = m:l
          Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
          Count = Count+1
          Peaklist[Count,] = NA
          if(Ratio<10|is.na(Gainloss_AIC)){
            Peaklist[(Count-1),] = NA
          }
          break
        } else{
          break
        }
      }
      level = quan_peak$level
      Peaklist$Start_RT[Count] = EIC$rt[m]
      Peaklist$End_RT[Count] = EIC$rt[l]
      Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
      Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
      Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
      Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
      if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
      } else {
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
      }
      Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                          error=function(e){"too_small"})
      Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
      Peaklist$RT[Count] = median(Peaktop$rt)
      Smoothie1 = Smoothie
      Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
      #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
      tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                               resp = Smoothie1$int,
                               bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
      
      RMSD_Ratio = Background/Baseline*100
      if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
        RMSD_Ratio = 0
      }
      if(RMSD_Ratio<10){
        Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
      } else {
        Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
      }
      
      if(is.na(Ratio)){
        Ratio = 0
      }
      Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
      
      if(grepl("GC",sample)){
        if(Differenci<=0.01&Ratio<3){
          Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
        }
      } else {
        if(Differenci<=0.05&Ratio<3){
          Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
        }
      }
      #Ratio of 1 % of Peak to Maximum has to be at least 10!
      
      #Constant = no real peak, too flat
      Constant_AIC = tcpl_fit$cnst_aic
      #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
      Gainloss_AIC = tcpl_fit$gnls_aic
      
      if(!is.na(Gainloss_AIC)){
        if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
          Gainloss_AIC = NA
        }
      }
      
      Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
      Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
      Peaklist$Height[Count] = max(EIC$int[m:l])
      id = m:l
      Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
      Start = l+1
      Count = Count+1
      Already = l
      Peaklist[Count,] = NA
      if(Ratio<10|is.na(Gainloss_AIC)){
        Peaklist[(Count-1),] = NA
      }
    }
    Peaklist = Peaklist[!is.na(Peaklist$RT)&Peaklist$Sequence!="too_small",]
    
    if(nrow(Peaklist)==0){
      Peaklist[1,] = NA
    }
    
    Check_sequence = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters)&Peaklist$Area > minimum_peak_area
    
    Filter = !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence
    Peaklist_final = filter(Peaklist,Filter==T)
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA, "LOD" = NA, "LOQ" = NA)
    }
  } else {
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA, "LOD" = NA, "LOQ" = NA)
  }
  Peaklist_final$Compound = quan_peak$Compound
  Peaklist_final$mz = quan_peak$mz
  
  if(!is.na(Peaklist_final$Height[1])){
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA, "LOD" = NA, "LOQ" = NA)
      Peaklist_final$Compound = quan_peak$Compound
      Peaklist_final$mz = quan_peak$mz
    }
  }
  
  if(length(Peaklist_final$RT[!is.na(Peaklist_final$RT)])>0){
    Peaklist_final = Peaklist_final[!is.na(Peaklist_final$RT),]
  }
  
  if(length(Peaklist_final$Comment)>1){
    ID = as.numeric(substr(Peaks$Compound[t],1,4))
    if(is.na(ID)){
      ID = paste0("IS",substr(Peaks$Compound[t],3,4))
    }
    if(ID%in%Multiple$ID){
      if(Multiple$multiple[match(ID,Multiple$ID)]==0){
        #first check if very small peak is included by mistake
        Peaklist_final = Peaklist_final[Peaklist_final$Height>inner_intensity_ratio_multiple_peaks*max(Peaklist_final$Height,na.rm=T),]
        
        potential_width = abs(Peaklist_final$End_RT_level[nrow(Peaklist_final)]-Peaklist_final$Start_RT_level[1])
        aimed_width = quan_peak$End_RT_level-quan_peak$Start_RT_level
        
        while(potential_width/aimed_width>1.2&nrow(Peaklist_final)>1){
          #Remove the peak which is most far away from center if widths are deviating (likely unspecific peak present)
          shift_corrected_times = Peaks$RT[t]+Peaklist_final$intensity_shift+IS_dependent_shift_compound
          center = shift_corrected_times[1]
          diffs = abs(Peaklist_final$RT-center)
          Peaklist_final = Peaklist_final[diffs!=max(diffs,na.rm = T),]
          potential_width = abs(Peaklist_final$End_RT_level[nrow(Peaklist_final)]-Peaklist_final$Start_RT_level[1])
        }
        sample_peak = Peaklist_final[1,]
        sample_peak$Start_RT = min(Peaklist_final$Start_RT,na.rm = T)
        sample_peak$RT = median(Peaklist_final$RT,na.rm = T)
        sample_peak$End_RT = max(Peaklist_final$End_RT,na.rm = T)
        sample_peak$Sequence = SAX_consensus(Peaklist_final$Sequence)
        sample_peak$Width = sample_peak$End_RT-sample_peak$Start_RT
        sample_peak$Height = max(Peaklist_final$Height,na.rm = T)
        sample_peak$Area = sum(Peaklist_final$Area,na.rm = T)
      } else {
        diff = abs(Peaklist_final$RT-Peaks$RT[t])
        Peaklist_final = Peaklist_final[match(min(diff),diff),]
        sample_peak = Peaklist_final
      }
    } else {
      diff = abs(Peaklist_final$RT-Peaks$RT[t])
      diff.grid = expand.grid(diff,diff)
      diff.grid$diff = abs(diff.grid$Var1-diff.grid$Var2)
      diff.grid = diff.grid[diff.grid$diff!=0,]
      differences = unique(diff.grid$diff)
      nearest_peaks = min(diff)==diff.grid$Var1[diff.grid$diff==min(diff.grid$diff)]
      if(min(differences)<0.05&T%in%nearest_peaks){
        minidist = rep(NA,nrow(Peaklist_final))
        for(mi in 1:length(minidist)){
          minidist[mi] = SAX_mindist(Peaklist_final$Sequence[mi],Peaks$Sequence[t],Peaklist_final$Nr_of_Points[mi],SAX_reference_table)
        }
        Peaklist_final = Peaklist_final[match(min(minidist),minidist),]
        sample_peak = Peaklist_final
      } else {
        Peaklist_final = Peaklist_final[match(min(diff),diff),]
        sample_peak = Peaklist_final
      }
    }
  } else {
    sample_peak = Peaklist_final
  }
  
  sample_peak$LOD = LOD
  sample_peak$LOQ = LOQ
  
  if(bad_chroma==F&nrow(sample_peak[!is.na(sample_peak$Height),])>0){
    for(qlk in 1:nrow(sample_peak)){
      Peak = EIC[EIC$rt>=sample_peak$Start_RT[qlk]&EIC$rt<=sample_peak$End_RT[qlk],]
    }
  }
  
  return(sample_peak)
}

#'sample_analysis
#'
#'Analysis of each sample. Takes shifts and confirming ions into account.
#'
#' @param sample_path path to the sample / file batch
#' @param results_path path to the output directory
#' @param sample the name of the sample / file batch
#' @param files names of the raw data files to be analyzed
#' @param Peaks results of peak annotation, data frame of peak characteristics from references
#' @param Peaks_list Peaks as list per reference
#' @param IS_dependent_shift a data frame linking samples and shifts of the internal standards
#' @param Intensity_dependent_shift a data frame listing compounds which show intensity-dependent shifts in calibrations
#' @param Multiple file defining peak patterns if multiple peaks occur in the m/z channel within RT_range
#' @param FullscanMS1 the scan filter used for analysis of the main peaks
#' @param precursor_window the m/z window of precursor masses included when searching for MS2 scans in .mzML files
#' @param inner_intensity_ratio_multiple_peaks the intensity ratio based on the maximum intensity within the inner RT range, used to filter valid peaks, usually a lower ratio
#' @param Scan_filters list of scan filters used for screening, mainly used for MS2 annotation
#' @param IS_Assignment data frame linking compounds to their assigned internal standards
#' @param compounds_shift_corrected the list of targets with retention times corrected
#' @param ppm_val allowed variability of target m/z in ppm
#' @param method_time the time of the MS measurement
#' @param is_centroidedrelevant for mzML files if scan is centroided or not
#' @param gen.plots boolean if plots should be generated
#' @param use.MINDIST boolean if MINDIST should be applied as peak selection criterion
#' @param SAX_reference_table the SAX distance matrix used for calculating MINDIST
#' @param alphabet_size the alphabet size used for SAX functions
#' @param zigzag_trigger_threshold defines the maximum fraction of the length of the EIC which can consist of zigzag patterns before triggering a more inclusive background calculation
#' @param normal_background_quantile the quantile used to define background intensities from the background EIC
#' @param higher_background_quantile the quantile used to define background intensities from the background EIC if the chromatogram is showing quality flags 
#' @param minimum_background_ratio the minimum ratio expected for a valid chromatogram of maximum intensity and the normal background quantile
#' @param extended_baseline_factor a ratio applied to the baseline to define values which are considered as baseline
#' @param width_smoothing a numeric value needed to define the density of smoothing for the EIC
#' @param width_factor_background the factor applied to the peak width to define a background window
#' @param sample_search_window the search window in minutes added to both ends of the expected retention time used to screen for the peak of interest
#' @param maximum_nr_of_a within SAX the letter a defines values with the lowest values and a valid peak shall only have a distinct number of these low values to avoid too flat peaks
#' @param minimum_nr_of_high_intensity_letters dependent on the selected SAX alphabet, the SAX sequence of a valid peak should have at least the respective number of letters defining high intensities
#' @param minimum_peak_area the value of a peak area which needs to be exceeded to be a valid peak
#' @param maximum_peak_width the maximum allowed peak width in minutes
#' @param minimum_nr_of_datapoints_per_peak the minimal number of data points needed to define a valid peak
#' @param maximum_allowed_shift the maximum shift allowed in minutes
#' @param maximum_allowed_shift_ratio the maximum fraction of IS-dependent shift of the allowed shift, only relevant if allowed shift exceeds maximum value
#' @param maximum_cutoff_intensity if the intensity of cutoffs which split peaks exceed this value, the cutoff is considered a peak-splitting cutoff and thereby enforced
#' @param minimum_confirming_peak_height the minimum height a confirming ion should exceed if expected
#' @param max_MINDIST a value used as reference to interpret if the MINDIST of two SAX sequences is significantly different
#' @param minimum_datapoints_per_sample_peak the minimal number of data points which should define a valid peak in target screening
#' @param t index of the compound within Peaks
#' @param j index of the file to be analyzed in files
#' @param rep boolean if the function was already called before
#' @export
sample_analysis = function(sample_path,results_path,sample,files,Peaks,Peaks_list,Solvent_Blank_ID,IS_dependent_shift,Intensity_dependent_shift,Multiple,FullscanMS1,precursor_window,inner_intensity_ratio_multiple_peaks,Scan_filters,IS_Assignment,compounds_shift_corrected,RT_range,ppm_val,method_time,IS_deviation,gen.plots,use.MINDIST,SAX_reference_table,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,minimum_background_ratio,extended_baseline_factor,width_smoothing,width_factor_background,sample_search_window,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,maximum_allowed_shift,maximum_allowed_shift_ratio,minimum_cutoff_intensity_factor,minimum_confirming_peak_height,max_MINDIST,minimum_datapoints_per_sample_peak,minimum_qualitative_threshold,t,j,rep){
  
  quan_peak = Peaks[t,]
  ID = as.numeric(substr(Peaks$Compound[t],1,4))
  if(is.na(ID)){
    ID = substr(Peaks$Compound[t],1,4)
  }
  
  Intensity_dependent_shift_names = names(Intensity_dependent_shift)
  
  if(compounds_shift_corrected$MS1.MS2[match(ID,compounds_shift_corrected$ID)]=="MS1"){
    MS_filter = FullscanMS1
  } else {
    diff = abs(Peaks$mz[t]-as.numeric(names(Scan_filters)))
    MS_filter = Scan_filters[[match(min(diff),diff)]]
  }
  
  if(is.null(quan_peak)){
    sample_peak = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
    sample_peak$Compound = Peaks$Compound[t]
    sample_peak$mz = Peaks$mz[t]
    
    if (T%in%grepl(".mzML", files)) {
      if(MS_filter==FullscanMS1){
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "EIC",
                               mz = sample_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
      } else {
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "MS2",
                               mz = sample_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        precursor = sample_peak$mz
        EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
        EIC = EIC[(EIC$mz>=RaMS::pmppm(sample_peak$mz,ppm_val)[1]&EIC$mz<=RaMS::pmppm(sample_peak$mz,ppm_val)[2]),]
        EIC = EIC[,c(1:2)]
      }
      
    } else {
      Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                       mass = sample_peak$mz,
                                       tol = ppm_val*2,
                                       filter = MS_filter,
                                       type = "xic")
      
      EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
    }
    
    ID_ = stringr::str_sub(Peaks$Compound[t],1,5)
    k = match(T,grepl(ID_,compounds_shift_corrected$identity))
    aimed_time = compounds_shift_corrected$`retention time`[k]
    test = EIC[EIC$rt>=(aimed_time-5*RT_range)&EIC$rt<=(aimed_time+5*RT_range),]
    if(nrow(test)>=5){
      EIC = test
    } else {
      if(length(EIC$rt)<=1){
        EIC = data.frame("rt"=c(0,method_time/3,method_time*2/3,method_time),"int"=c(0,0,0,0))
      }
    }
    if(gen.plots==T){
      dir.create(paste0(results_path,"/plots"))
      if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
        Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
        PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      } else {
        PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      }
      pdf(file=PDF_path)
      p <- ggplot(EIC, aes(rt, int)) + 
        geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
      p = p + 
        labs(title = stringr::str_wrap(paste0(files[j]," - ",Peaks$Compound[t]),width=70), 
             subtitle = "NO_VALID_REFERENCE") +
        theme_classic() +
        theme(plot.title = element_text(lineheight = 1.1),
              plot.subtitle = element_text(lineheight = 1.0,colour = "darkgrey"),
              legend.position = "none")
      p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                               labels = scales::number_format(accuracy = 0.01))
      print(p)
      dev.off()
    }
    return(sample_peak)
  }
  
  if(is.nan(quan_peak$RT)|is.na(quan_peak$RT)){
    sample_peak = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
    sample_peak$Compound = Peaks$Compound[t]
    sample_peak$mz = Peaks$mz[t]
    if (T%in%grepl(".mzML", files)) {
      if(MS_filter==FullscanMS1){
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "EIC",
                               mz = sample_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
      } else {
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "MS2",
                               mz = sample_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        precursor = sample_peak$mz
        EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
        EIC = EIC[(EIC$mz>=RaMS::pmppm(sample_peak$mz,ppm_val)[1]&EIC$mz<=RaMS::pmppm(sample_peak$mz,ppm_val)[2]),]
        EIC = EIC[,c(1:2)]
      }
    } else {
      
      Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                       mass = sample_peak$mz,
                                       tol = ppm_val*2,
                                       filter = MS_filter,
                                       type = "xic")
      
      EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
    }
    ID_ = stringr::str_sub(Peaks$Compound[t],1,5)
    k = match(T,grepl(ID_,compounds_shift_corrected$identity))
    aimed_time = compounds_shift_corrected$`retention time`[k]
    test = EIC[EIC$rt>=(aimed_time-5*RT_range)&EIC$rt<=(aimed_time+5*RT_range),]
    if(nrow(test)>=5){
      EIC = test
    } else {
      if(length(EIC$rt)<=1){
        EIC = data.frame("rt"=c(0,method_time/3,method_time*2/3,method_time),"int"=c(0,0,0,0))
      }
    }
    if(gen.plots==T){
      dir.create(paste0(results_path,"/plots"))
      if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
        Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
        PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      } else {
        PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      }
      pdf(file=PDF_path)
      p <- ggplot(EIC, aes(rt, int)) + 
        geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
      p = p + 
        labs(title = stringr::str_wrap(paste0(files[j]," - ",Peaks$Compound[t]),width=70), 
             subtitle = "NOT_VALID_REFERENCE") +
        theme_classic() +
        theme(plot.title = element_text(lineheight = 1.1),
              plot.subtitle = element_text(lineheight = 1.0,colour = "darkgrey"),
              legend.position = "none")
      p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                               labels = scales::number_format(accuracy = 0.01))
      print(p)
      dev.off()
    }
    return(sample_peak)
  }
  
  if (T%in%grepl(".mzML", files)) {
    if(MS_filter==FullscanMS1){
      EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                             grab_what = "EIC",
                             mz = quan_peak$mz,
                             ppm = ppm_val,
                             rtrange = NULL,
                             prefilter = -1)
      EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
    } else {
      EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                             grab_what = "MS2",
                             mz = quan_peak$mz,
                             ppm = ppm_val,
                             rtrange = NULL,
                             prefilter = -1)
      precursor = quan_peak$mz
      EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
      EIC = EIC[(EIC$mz>=RaMS::pmppm(quan_peak$mz,ppm_val)[1]&EIC$mz<=RaMS::pmppm(quan_peak$mz,ppm_val)[2]),]
      EIC = EIC[,c(1:2)]
    }
  } else {
    
    Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                     mass = quan_peak$mz,
                                     tol = ppm_val*2,
                                     filter = MS_filter,
                                     type = "xic")
    
    EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
  }
  
  EIC$int[EIC$int==Inf] = 0
  
  #Calculate shift based on intensity and IS-shift
  IS_Assignment_shortened = IS_Assignment
  IS_Assignment_shortened$Compound = sapply(IS_Assignment$Compound,Remove_ion)
  IS_Assignment_shortened$ISTD = sapply(IS_Assignment$ISTD,Remove_ion)
  if(is.na(as.numeric(ID))){
    ISTD = IS_Assignment_shortened$ISTD[match(ID,sapply(IS_Assignment_shortened$Compound,function(x){stringr::str_sub(x,1,4)}))]
  } else {
    ISTD = IS_Assignment_shortened$ISTD[match(ID,sapply(IS_Assignment_shortened$Compound,function(x){as.numeric(stringr::str_sub(x,1,4))}))]
  }
  
  if(is.na(ISTD)){
    ISTD = Peaks$Compound[t]
  }
  IS_dependent_shift_compound = IS_dependent_shift[match(ISTD,IS_dependent_shift$ISTD),(j+1)]
  potential_intensity_dependent_shift = 0
  if(Peaks$Compound[t]%in%Intensity_dependent_shift_names){
    potential_intensity_dependent_shift = max(abs(Intensity_dependent_shift[[match(Peaks$Compound[t],Intensity_dependent_shift_names)]]$shift),na.rm = T)
  }
  allowed_shift = abs(IS_dependent_shift_compound)+potential_intensity_dependent_shift
  
  if(allowed_shift>0.5){
    if(T%in%grepl("pos|neg",sample)){
      width_smoothing = 3/0.5*allowed_shift
    } else {
      width_smoothing = 3/0.5*allowed_shift*10
    }
  } else {
    width_smoothing = 3 #empirical value
  }
  
  Background_window = quan_peak$Width*width_factor_background
  
  if(Background_window<width_factor_background*sample_search_window){
    Background_window = width_factor_background*sample_search_window
  }
  
  Filter_background =  EIC$rt > (Peaks$Start_RT[t]-Background_window-allowed_shift) & EIC$rt < (Peaks$End_RT[t]+Background_window+allowed_shift)
  
  width = abs(quan_peak$End_RT_level-quan_peak$Start_RT_level)
  
  Filter = EIC$rt < (Peaks$Start_RT_level[t]-sample_search_window-width-allowed_shift) | EIC$rt > (Peaks$End_RT_level[t]+sample_search_window+width+allowed_shift)
  
  bad_chroma = F
  
  EIC_background = dplyr::filter(EIC,Filter_background==T)
  EIC = dplyr::filter(EIC,Filter==F)
  if(length(EIC$rt)<3|length(EIC_background$rt)<3){
    bad_chroma = T
  }
  if(bad_chroma==F){
    plot(EIC$int~EIC$rt,type="p",main=paste0(quan_peak$Compound))
  }
  timesplit = ceiling(length(EIC$rt)/width_smoothing)
  if(bad_chroma==F){
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(i in 1:width_smoothing){
      EICs[[i]] = EIC[(Start:(Start+timesplit)),]
      EICs[[i]] = EICs[[i]][!is.na(EICs[[i]]$rt)&!is.na(EICs[[i]]$int),]
      timerange = max(EICs[[i]]$rt)-min(EICs[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_smooth[[i]] = as.data.frame(ksmooth(EICs[[i]]$rt,EICs[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[i]] = EICs_smooth[[i]][!is.na(EICs_smooth[[i]]$x)&!is.na(EICs_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    lines(EIC_smooth$rt,EIC_smooth$int,col="red",type = "p")
    
    threshold = as.numeric(quantile(EIC_background$int,probs=c(normal_background_quantile)))
    for(we in 1:nrow(EIC_background)){
      if(we == 1){
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = "b"
        } else {
          sequenci = "a"
        }
      } else {
        if(!is.na(threshold)&EIC_background$int[we]>threshold){
          sequenci = append(sequenci,"b")
        } else {
          sequenci = append(sequenci,"a")
        }
      }
    }
    sequencis = paste0(sequenci,collapse = "")
    
    EICs_background = list()
    EICs_background_smooth = list()
    timesplit = ceiling(length(EIC_background$rt)/width_smoothing)
    Start = 1
    for(i in 1:width_smoothing){
      EICs_background[[i]] = EIC_background[(Start:(Start+timesplit)),]
      EICs_background[[i]] = EICs_background[[i]][!is.na(EICs_background[[i]]$rt)&!is.na(EICs_background[[i]]$int),]
      timerange = max(EICs_background[[i]]$rt)-min(EICs_background[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_background_smooth[[i]] = as.data.frame(ksmooth(EICs_background[[i]]$rt,EICs_background[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_background_smooth[[i]] = EICs_background_smooth[[i]][!is.na(EICs_background_smooth[[i]]$x)&!is.na(EICs_background_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC_background$rt)){
        break
      }
    }
    
    EIC_background_smooth = dplyr::bind_rows(EICs_background_smooth)
    colnames(EIC_background_smooth) = c("rt","int")
    zigzag_positions = as.data.frame(stringr::str_locate_all(sequencis,"aba"))
    zigzag_positions = rbind(zigzag_positions,as.data.frame(stringr::str_locate_all(sequencis,"abba")))
    zigzag_positions = rbind(zigzag_positions,as.data.frame(stringr::str_locate_all(sequencis,"abbba")))
    
    if(nrow(zigzag_positions)>0){
      for(zi in 1:nrow(zigzag_positions)){
        Start = zigzag_positions$start[zi]
        End = zigzag_positions$end[zi]
        sequenci[Start:End] = "c"
      }
    }
    
    abline(v=EIC_background$rt[sequenci=="c"],col=adjustcolor(col="lightgrey",alpha.f = 0.2))
    
    if (T%in%grepl(".mzML", files)){
      mean_Background = as.numeric(quantile(EIC_background$int,na.rm = T,probs = c(normal_background_quantile)))
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    } else {
      mean_Background = mean(EIC_background$int,na.rm = T)
      max_Background = max(EIC_background$int,na.rm = T)
      if(length(mean_Background)>0&length(max_Background)>0&mean_Background!=0){
        Background_ratio = max_Background/mean_Background
      } else {
        Background_ratio = 0
      }
    }
    
    if((length(sequenci[sequenci=="c"])>(zigzag_trigger_threshold*length(sequenci)))&(max(EIC_background$int,na.rm = T)!=min(EIC_background$int,na.rm=T))){
      Q = as.numeric(quantile(EIC_background$int[sequenci=="c"],probs=c(higher_background_quantile)))
      Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
      Background = sd(Background_signals,na.rm=T)
      Baseline = mean(Background_signals)
    } else if(nrow(EIC_smooth[EIC_smooth$int!=0,])>0){
      if(Background_ratio<minimum_background_ratio){
        Q = as.numeric(quantile(EIC_background$int,probs=c(higher_background_quantile)))
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<Q]
          Background = sd(Background_signals,na.rm=T)
          Baseline = mean(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
        
      } else {
        Q = as.numeric(quantile(EIC_background_smooth$int,probs=c(normal_background_quantile)))
        Q = ifelse(is.na(Q),0,Q)
        if(Q!=0){
          Background_signals = EIC_background_smooth$int[EIC_background_smooth$int<=Q]
          Baseline = mean(Background_signals,na.rm=T)
          Background = sd(Background_signals,na.rm=T)
        } else {
          Background = 0
          Baseline = 0
        }
      }
    } else {
      Background = 0
      Baseline = 0
    }
    
    nrow(EIC_smooth[EIC_smooth$int<=extended_baseline_factor*min(EIC_smooth$int),])>(zigzag_trigger_threshold*length(sequenci))
    
    abline(h=Baseline,col="blue",lwd=3)
    
  } else {
    EIC_smooth = EIC
    EIC_background = data.frame("rt"=c(1,2,3),"int"=c(0.000001,0.000001,0.000001))
    Background = 0
    Baseline = 0
  }
  
  if(length(EIC$rt)<=1){
    EIC = data.frame("rt"=c(1,2,3),"int"=c(0,0,0))
  }
  
  LOD = as.numeric(Baseline+3*Background)
  LOQ = as.numeric(Baseline+10*Background)
  
  if(is.na(LOD)|is.nan(LOD)){
    bad_chroma = T
  }
  
  if(bad_chroma==F){
    abline(h=LOD,col="green")
    abline(h=LOQ,col="orange")
  }
  guggl = seewave::SAX(EIC$int,alphabet_size = 3,PAA_number = 10)
  
  if(!(T%in%grepl("c",guggl))){
    bad_chroma = F
  }
  
  diffl = abs(EIC_smooth$rt-(min(EIC_smooth$rt,na.rm = T)+0.1))
  density = round(match(min(diffl,na.rm = T),diffl),digits = 0)
  
  SAX_result = rep("a",length(EIC$int))
  
  SAX_result_smooth = rep("a",length(EIC_smooth$int))
  
  SAX_result_small = SAX_result
  
  SAX_result_shift = SAX_result
  
  if(bad_chroma==F){
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>LOD){
        SAX_result_small[l] = "b"
      } else {
        SAX_result_small[l] = "a"
      }
    }
    
    for(l in 1:length(EIC$int)){
      if(EIC$int[l]>(Baseline*extended_baseline_factor)){
        SAX_result_shift[l] = "b"
      } else {
        SAX_result_shift[l] = "a"
      }
    }
    
    if(length(EIC$int)<3|var(EIC$int)==0){
      bad_chroma = T
    }
    
    if(!(T%in%grepl("b",SAX_result_small))){
      bad_chroma = T
    }
  }
  
  if(bad_chroma == F){
    Peaklist = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA)
    
    if(T%in%grepl("b",SAX_result_small)){
      Peak = match(max(EIC_smooth$int,na.rm = T),EIC_smooth$int)
      Peak_rt = EIC_smooth$rt[Peak]
      Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
      Cutoffs = c()
      if(!is.na(Peak)&Peak_int>LOD){
        remaining_Peaks = T
      } else {
        remaining_Peaks = F
      }
      while(remaining_Peaks == T){
        Peak_start = NA
        Peak_end = NA
        Cutoff_r = NA
        Cutoff_l = NA
        if(Peak>1){
          Begin_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak-1):1]
          Points_index = (Peak-1):1
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_start = Points_index[c-counter]
                break
              }
              Peak_start = Points_index[[c]]+1
              break
            }
            if(Check[c]>=(Begin_int*0.7)){
              SAX_result_smooth[Points_index[c]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_start = Points_index[c-7]
              if(c==7){
                Peak_start = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_start = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_start = 1
                } else {
                  Peak_start = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_start = 1
        }
        
        if(Peak < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak]
          Points_int = EIC_smooth$int[(Peak+1):length(SAX_result_smooth)]
          Points_index = (Peak+1):length(SAX_result_smooth)
          Check = Points_int
          counter = 0
          for(c in 1:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Peak_end = Points_index[c-counter]
                break
              }
              Peak_end = Points_index[[c]]-1
              break
            }
            if(Check[c]>=(0.7*End_int)){
              SAX_result_smooth[Points_index[c]] = "b"
              SAX_result_smooth[Points_index[c-1]] = "b"
              while(counter!=0){
                SAX_result_smooth[Points_index[c-counter]] = "b"
                counter = counter-1
              }
            } else {
              counter = counter+1
            }
            if(counter==7){
              Peak_end = Points_index[c-7]
              if(c==7){
                Peak_end = Peak
              }
              break
            }
            if(c==length(Points_index)){
              if(counter==0){
                Peak_end = Points_index[length(Points_index)]
              } else {
                if(c==counter){
                  Peak_end = length(SAX_result_smooth)
                } else {
                  Peak_end = Points_index[length(Points_index)-counter]
                }
              }
            }
          }
        } else {
          Peak_end = length(SAX_result_smooth)
        }
        
        SAX_result_smooth[Peak] = "b"
        
        if(Peak_start > 1){
          Begin_int = EIC_smooth$int[Peak_start]
          Points_int = EIC_smooth$int[(Peak_start):1]
          Points_index = (Peak_start):1
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_l = Points_index[c-counter]
                break
              }
              Cutoff_l = Points_index[c-1]
              break
            }
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_l = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_l = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              if(counter==(5*density+1)){
                Cutoff_l = Points_index[(c-((5*density+1)-1))]
                break
              }
            }
          }
          if(is.na(Cutoff_l)){
            Cutoff_l = 1
          }
        } else {
          Cutoff_l = 1
        }
        
        if(Peak_end < length(SAX_result_smooth)){
          End_int = EIC_smooth$int[Peak_end]
          Points_int = EIC_smooth$int[Peak_end:length(SAX_result_smooth)]
          Points_index = Peak_end:length(SAX_result_smooth)
          Check = Points_int
          counter = 0 #counts negative trends
          counter2 = 0 #counts positive trends
          threshold = NA
          for(c in 2:length(Check)){
            if(SAX_result_smooth[Points_index[c]]=="b"){
              if(counter!=0){
                Cutoff_r = Points_index[c-counter]
                break
              }
              Cutoff_r = Points_index[c-1]
              break
            }
            if(counter==0){
              if(Check[c]<(Check[c-1])){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
              } else {
                counter = counter+1
                if(is.na(threshold)){
                  threshold = Points_int[c-1]*0.6
                }
                Cutoff_r = Points_index[c-1]
              }
            }else{
              if(Check[c]>(1.4*1.6666666*threshold)){
                counter2 = counter2+1
              }
              if(counter2>(0.25*density+1)){
                threshold=0
              }
              if(Check[c]<threshold){
                SAX_result_smooth[Points_index[c]] = "b"
                SAX_result_smooth[Points_index[c-1]] = "b"
                Cutoff_r = Points_index[c]
                threshold = NA
                counter2 = 0
                while(counter!=0){
                  SAX_result_smooth[Points_index[c-counter]] = "b"
                  counter = counter-1
                }
              } else {
                counter = counter+1
              }
              
              if(counter==(5*density+1)){
                Cutoff_r = Points_index[(c-((5*density+1)-1))]
                break
              }
            }
          }
          if(is.na(Cutoff_r)){
            Cutoff_r = length(EIC_smooth$rt)
          }
        } else {
          Cutoff_r = length(EIC_smooth$rt)
        }
        
        Pairs = length(Cutoffs)/2
        if(Pairs>0){
          for(pa in 1:Pairs){
            pairis = data.frame("Start"=seq(1,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2),"End"=seq(2,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2))
            Range = Cutoffs[pairis$Start[pa]]:Cutoffs[pairis$End[pa]]
            if(Cutoff_l%in%Range){
              Cutoff_l = max(Range)
            }
            if(Cutoff_r%in%Range){
              Cutoff_r = min(Range)
            }
          }
        }
        
        Cutoffs = append(Cutoffs,c(Cutoff_l,Cutoff_r))
        
        Indizes = 1:length(EIC_smooth$int)
        
        Filter = SAX_result_smooth!="b"
        
        Indizes = Indizes[Filter]
        
        Peak = Indizes[match(max(EIC_smooth$int[SAX_result_smooth!="b"]),EIC_smooth$int[Filter])]
        
        Peak_rt = EIC_smooth$rt[Peak]
        
        Peak_int = max(EIC$int[EIC$rt>(Peak_rt-0.05)&EIC$rt<(Peak_rt+0.05)])
        
        if(!is.na(Peak)&Peak_int>LOD){
          remaining_Peaks = T
        } else {
          remaining_Peaks = F
        }
      }
      SAX_result_smooth[Cutoffs] = "a"
      if(1%in%Cutoffs){
        SAX_result_smooth[1] = "b"
      }
      if(length(SAX_result_smooth)%in%Cutoffs){
        SAX_result_smooth[length(SAX_result_smooth)] = "b"
      }
    }
    
    Start = 1
    
    while(Start<length(SAX_result_smooth)){
      peakbegin = match("b",SAX_result_smooth[Start:length(SAX_result_smooth)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result_smooth[peakbegin:length(SAX_result_smooth)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result_smooth)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      Starttime = EIC_smooth$rt[peakbegin]
      Endtime = EIC_smooth$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    SAX_result[EIC$int<=extended_baseline_factor*Baseline] = "a"
    
    Start = 1
    while(Start<length(SAX_result)){
      peakbegin = match("b",SAX_result[Start:length(SAX_result)])+Start-1
      if(is.na(peakbegin)){
        break
      }
      peakend = match("a",SAX_result[peakbegin:length(SAX_result)])+peakbegin-2
      if(is.na(peakend)){
        peakend = length(SAX_result)
        if(peakbegin==peakend){
          break
        }
      }
      if(peakbegin==peakend){
        Start = peakend+1
        next
      }
      if(peakbegin>1){
        if(EIC$int[peakbegin-1]<=(extended_baseline_factor*Baseline)){
          peakbegin = peakbegin - 1
        }
      }
      if(peakend<length(SAX_result)){
        if(EIC$int[peakend+1]<=(extended_baseline_factor*Baseline)){
          peakend = peakend + 1
        }
      }
      Starttime = EIC$rt[peakbegin]
      Endtime = EIC$rt[peakend]
      SAX_result[EIC$rt>=Starttime&EIC$rt<=Endtime] = "b"
      Start = peakend + 1
    }
    
    Cutoffs_time = EIC_smooth$rt[Cutoffs]
    
    for(klo in 1:length(Cutoffs_time)){
      diff = abs(EIC$rt-Cutoffs_time[klo])
      time = EIC$rt[match(min(diff,na.rm = T),diff)]
      SAX_result[EIC$rt==time] = "a"
    }
    
    #Test whether it is plausible that there are traces of analyte or just non-specific nearby peaks
    if(is.null(Cutoffs)){
      Cutoffs = NA
    }
    if(!is.na(Cutoffs)&length(Cutoffs)>0){
      n_maxima = length(Cutoffs)/2
      cutoff_maxima = rep(NA,n_maxima)
      cutoff_indices = matrix(c(1:length(Cutoffs)),ncol=2,byrow = T)
      for(c in 1:n_maxima){
        cutoff_maxima[c] = max(EIC_smooth$int[Cutoffs[cutoff_indices[c,1]]:Cutoffs[cutoff_indices[c,2]]])
      }
      cutoff_maxima_rt = EIC_smooth$rt[which(EIC_smooth$int%in%cutoff_maxima)]
      integrate_peak_rt = F
      for(d in 1:n_maxima){
        if(allowed_shift>quan_peak$RT_tol){
          cutoff_rt = ifelse(cutoff_maxima_rt[d]<=quan_peak$RT+allowed_shift&cutoff_maxima_rt[d]>=quan_peak$RT-allowed_shift,cutoff_maxima_rt[d],NA)
        } else {
          cutoff_rt = ifelse(cutoff_maxima_rt[d]<=quan_peak$RT+quan_peak$RT_tol&cutoff_maxima_rt[d]>=quan_peak$RT-quan_peak$RT_tol,cutoff_maxima_rt[d],NA)
        }
        if(is.na(cutoff_rt)){
          next
        }
        diff_cutoff = abs(cutoff_rt-cutoff_maxima_rt)
        if((min(diff_cutoff[diff_cutoff!=0])<quan_peak$Width*width_factor_background)){
          integrate_peak_rt = T
          break
        } 
      }
    } else {
      integrate_peak_rt = F
    }
    
    if(integrate_peak_rt){
      if(quan_peak$Start_RT_level==quan_peak$End_RT_level){
        EIC_sub = EIC[EIC$rt>=quan_peak$Start_RT&EIC$rt<=quan_peak$End_RT,]
        EIC_sub_narrow = EIC[EIC$rt>=quan_peak$Start_RT_level&EIC$rt<=quan_peak$End_RT_level,]
        EIC_not_peak = EIC_smooth[EIC_smooth$rt%in%cutoff_maxima_rt[!(cutoff_maxima_rt>=min(EIC_sub$rt)&cutoff_maxima_rt<=max(EIC_sub$rt))],]
        if(nrow(EIC_not_peak)>0){
          cutoff_background_ratio = quantile(EIC_sub$int,probs = c(0.95),na.rm = T)/quantile(EIC_not_peak$int,probs = c(0.95),na.rm = T)
          if(is.na(cutoff_background_ratio)|is.nan(cutoff_background_ratio)){
            cutoff_background_ratio = 0
          }
        } else {
          cutoff_background_ratio = minimum_background_ratio
        }
        intensities_in_centre = nrow(EIC_sub_narrow[EIC_sub_narrow$int>LOD,])>0
        if((cutoff_background_ratio<minimum_background_ratio)&intensities_in_centre){
          SAX_aimed_peak = rep("a",nrow(EIC_sub))
        } else {
          SAX_aimed_peak = SAX_result[EIC$rt>=quan_peak$Start_RT&EIC$rt<=quan_peak$End_RT]
        }
      } else {
        EIC_sub = EIC[EIC$rt>=quan_peak$Start_RT&EIC$rt<=quan_peak$End_RT,]
        EIC_sub_narrow = EIC[EIC$rt>=quan_peak$Start_RT_level&EIC$rt<=quan_peak$End_RT_level,]
        EIC_not_peak = EIC_smooth[EIC_smooth$rt%in%cutoff_maxima_rt[!(cutoff_maxima_rt>=min(EIC_sub$rt)&cutoff_maxima_rt<=max(EIC_sub$rt))],]
        if(nrow(EIC_not_peak)>0){
          cutoff_background_ratio = quantile(EIC_sub$int,probs = c(0.95),na.rm = T)/quantile(EIC_not_peak$int,probs = c(0.95),na.rm = T)
          if(is.na(cutoff_background_ratio)|is.nan(cutoff_background_ratio)){
            cutoff_background_ratio = 0
          }
        } else {
          cutoff_background_ratio = minimum_background_ratio
        }
        intensities_in_centre = nrow(EIC_sub_narrow[EIC_sub_narrow$int>LOD,])>0
        if((cutoff_background_ratio<minimum_background_ratio)&intensities_in_centre){
          SAX_aimed_peak = rep("a",nrow(EIC_sub))
        } else {
          SAX_aimed_peak = SAX_result[EIC$rt>=quan_peak$Start_RT&EIC$rt<=quan_peak$End_RT]
        }
      }
      
      if(T%in%grepl("b",SAX_aimed_peak)){
        Start = match("b",SAX_aimed_peak)
        End = length(SAX_aimed_peak)-match("b",rev(SAX_aimed_peak))+1
        SAX_aimed_peak[Start:End] = "b"
        SAX_result[EIC$rt>=quan_peak$Start_RT&EIC$rt<=quan_peak$End_RT] = SAX_aimed_peak
      }
    }
    
    #If Cutoffs are separating two peaks which are close to each other, apply them if > LOQ
    Cutoffs_time = EIC_smooth$rt[Cutoffs]
    Cutoffs_intensity = EIC_smooth$int[Cutoffs]
    if(LOQ>0){
      Cutoffs_time = Cutoffs_time[Cutoffs_intensity>LOQ]
    } else {
      Cutoffs_time = Cutoffs_time[Cutoffs_intensity>minimum_cutoff_intensity_factor*quan_peak$Height]
    }
    
    for(klo in 1:length(Cutoffs_time)){
      diff = abs(EIC$rt-Cutoffs_time[klo])
      time = EIC$rt[match(min(diff,na.rm = T),diff)]
      SAX_result[EIC$rt==time] = "a"
    }
    
    Start = 1
    End = length(SAX_result)
    Count = 1
    Already = 0
    
    while(5==5){
      not_a = match(letters[2:26],SAX_result[Start:End])+Already
      k = not_a[order(not_a,decreasing = F)[1]]
      k = k 
      if(!is.na(k)){
        k = k
      }
      l = match("a",SAX_result[k+1:End])+k-Start+Already+1
      if(!is.na(l)){
        l = l
      }
      #test for gaps
      m = k - 1
      if(is.na(m)){
        break
      }
      if(m==0){
        m = 1
      }
      if(abs(EIC$rt[k]-EIC$rt[m])>0.01){
        m = k
      }
      if(is.na(l)){
        if(mean(!is.na(not_a))!=0){ 
          level = quan_peak$level
          l = End
          Peaklist$Start_RT[Count] = EIC$rt[m]
          Peaklist$End_RT[Count] = EIC$rt[l]
          Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
          Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
          Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
          Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
          if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
          } else {
            Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
          }
          Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                              error=function(e){"too_small"})
          Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
          Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
          Peaklist$Height[Count] = max(EIC$int[m:l])
          
          Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
          Peaklist$RT[Count] = median(Peaktop$rt)
          
          Smoothie1 = Smoothie
          Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
          #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
          tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                                   resp = Smoothie1$int,
                                   bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
          
          RMSD_Ratio = Background/Baseline*100
          if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
            RMSD_Ratio = 0
          }
          if(RMSD_Ratio<10){
            Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
          } else {
            Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
          }
          
          if(is.na(Ratio)){
            Ratio = 0
          }
          Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
          
          Threshold_Differenci = plyr::round_any(Differenci,0.01,f=floor)
          Threshold_Ratio = plyr::round_any(Ratio,0.1,f=ceiling)
          Threshold_Accepted = Threshold_Differenci<=0.1&Threshold_Ratio>=2
          if(grepl("GC",sample)){
            if((Differenci<=0.01&Ratio<3)|Threshold_Accepted){
              Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
            }
          } else {
            if((Differenci<=0.05&Ratio<3)|Threshold_Accepted){
              Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
            }
          }
          
          #Constant = no real peak, too flat
          Constant_AIC = tcpl_fit$cnst_aic
          #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
          Gainloss_AIC = tcpl_fit$gnls_aic
          
          if(!is.na(Gainloss_AIC)){
            if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
              Gainloss_AIC = NA
            }
          }
          #If gnls_aic = NA or is bigger than Constant_AIC, there is likely NO valid peak
          #Peaklist$RT[Count] = EIC_smooth$rt[match(max(EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]]),EIC_smooth$int[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l]])+match(T,EIC_smooth$rt>=EIC$rt[m])-1]
          id = m:l
          Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
          Count = Count+1
          Peaklist[Count,] = NA
          if(Ratio<10|is.na(Gainloss_AIC)){
            Peaklist[(Count-1),] = NA
          }
          break
        } else{
          break
        }
      }
      level = quan_peak$level
      Peaklist$Start_RT[Count] = EIC$rt[m]
      Peaklist$End_RT[Count] = EIC$rt[l]
      Smoothie = EIC_smooth[EIC_smooth$rt>=EIC$rt[m]&EIC_smooth$rt<=EIC$rt[l],]
      Peakl = EIC[EIC$rt>=EIC$rt[m]&EIC$rt<=EIC$rt[l],]
      Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
      Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
      if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
      } else {
        Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
      }
      Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                          error=function(e){"too_small"})
      Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
      Peaklist$RT[Count] = median(Peaktop$rt)
      
      Smoothie1 = Smoothie
      Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
      #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
      tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                               resp = Smoothie1$int,
                               bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
      
      RMSD_Ratio = plyr::round_any(Background/Baseline*100,1,ceiling)
      if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
        RMSD_Ratio = 0
      }
      if(RMSD_Ratio<10){
        Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
      } else {
        Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
      }
      
      if(is.na(Ratio)){
        Ratio = 0
      }
      Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
      Threshold_Differenci = plyr::round_any(Differenci,0.01,f=floor)
      Threshold_Ratio = plyr::round_any(Ratio,0.1,f=ceiling)
      Threshold_Accepted = Threshold_Differenci<=0.1&Threshold_Ratio>=2
      if(grepl("GC",sample)){
        if((Differenci<=0.01&Ratio<3)|Threshold_Accepted){
          Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
        }
      } else {
        if((Differenci<=0.05&Ratio<3)|Threshold_Accepted){
          Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
        }
      }
      
      #Constant = no real peak, too flat
      Constant_AIC = tcpl_fit$cnst_aic
      #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
      Gainloss_AIC = tcpl_fit$gnls_aic
      
      if(!is.na(Gainloss_AIC)){
        if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
          Gainloss_AIC = NA
        }
      }
      #If gnls_aic = NA or is bigger than Constant_AIC, there is likely NO valid peak
      Peaklist$Nr_of_Points[Count] = length(EIC$int[m:l])
      Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
      Peaklist$Height[Count] = max(EIC$int[m:l])
      id = m:l
      Peaklist$Area[Count] = sum(diff(EIC$rt[id])*rollmean(EIC$int[id],2))
      Start = l+1
      Count = Count+1
      Already = l
      Peaklist[Count,] = NA
      if(Ratio<3|is.na(Gainloss_AIC)){
        Peaklist[(Count-1),] = NA
      }
    }
    Peaklist = Peaklist[!is.na(Peaklist$RT)&Peaklist$Sequence!="too_small",]
    
    if(ID%in%Multiple$ID&nrow(Peaklist)>0){
      if(Multiple$multiple[match(ID,Multiple$ID)]==0){
        for(m in nrow(Peaklist)){
          Peaklist$RT[m] = median(c(Peaklist$Start_RT_level[m],Peaklist$End_RT_level[m]),na.rm = T)
        }
      }
    }
    
    if(nrow(Peaklist)==0){
      Peaklist[1,] = NA
    }
    
    Check_sequence = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters|Peaklist$Nr_of_Points>minimum_nr_of_datapoints_per_peak)&Peaklist$Area > minimum_peak_area
    
    Filter = !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence
    Peaklist_final = filter(Peaklist,Filter==T)
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
    }
  } else {
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
  }
  Peaklist_final$Compound = quan_peak$Compound
  Peaklist_final$mz = quan_peak$mz
  Peaklist_final$intensity_shift = 0
  
  if(!is.na(Peaklist_final$Height[1])){
    if(mean(Filter)==0 | length(Filter)==0){
      Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
      Peaklist_final$Compound = quan_peak$Compound
      Peaklist_final$mz = quan_peak$mz
    }
  }
  
  if(length(Peaklist_final$RT[!is.na(Peaklist_final$RT)])>0){
    Peaklist_final = Peaklist_final[!is.na(Peaklist_final$RT),]
    if(Peaks$Compound[t]%in%Intensity_dependent_shift_names){
      for(tz in 1:nrow(Peaklist_final)){
        frame_intensity_shift = Intensity_dependent_shift[[match(Peaks$Compound[t],Intensity_dependent_shift_names)]]
        frame_intensity_shift$diff = abs(Peaklist_final$Height[tz]-frame_intensity_shift$int)
        frame_intensity_shift$diff[frame_intensity_shift$diff==0] = 10^-100
        frame_intensity_shift$diff = 1/frame_intensity_shift$diff
        two_smallest = frame_intensity_shift[order(frame_intensity_shift$diff,decreasing = T),]
        total = sum(two_smallest$diff)
        intensity_dependent_shift = sum(two_smallest$shift*two_smallest$diff)/total
        Peaklist_final$intensity_shift[tz] = intensity_dependent_shift
      }
    }
  }
  #If too many peaks are found, use only the 10 nearest (takes too long for too many peaks with qual check)
  if(nrow(Peaklist_final[!is.na(Peaklist_final$RT),])>10){
    diff = abs(Peaklist_final$RT-(quan_peak$RT+allowed_shift))
    test = Peaklist_final[order(diff,decreasing = F)[1:10],]
    if(allowed_shift>maximum_allowed_shift&(IS_dependent_shift_compound/allowed_shift)>maximum_allowed_shift_ratio){
      time_max = Peaklist_final$RT[match(max(Peaklist_final$Area),Peaklist_final$Area)]
      if(time_max>=(quan_peak$RT-RT_range)&time_max<=(quan_peak$RT+RT_range)){
        #Assumption that allowed shift is wrong
        IS_dependent_shift_compound = 0
        allowed_shift = 0+potential_intensity_dependent_shift
        diff = abs(Peaklist_final$RT-(quan_peak$RT+allowed_shift))
        test = Peaklist_final[order(diff,decreasing = F)[1:10],]
      }
    }
    Peaklist_final = test
    Peaklist_final = Peaklist_final[!is.na(Peaklist_final$RT),]
  }
  
  Peaklist_final$confirmed = 0
  if(nrow(Peaklist_final[!is.na(Peaklist_final$RT),])>=1){
    for(quani in 1:nrow(Peaklist_final)){
      sample_peak = Peaklist_final[quani,]
      if(length(sample_peak$RT[!is.na(sample_peak$RT)])>0){
        #########################################################
        #Search for confirming ions - but only ones with expected
        #########################################################
        qual_peaks = Peaks_list[[t]][Peaks_list[[t]]$type!="quan",]
        if(length(qual_peaks$Compound)>0){
          for(q in 1:length(qual_peaks$Compound)){
            if(compounds_shift_corrected$MS1.MS2[match(ID,compounds_shift_corrected$ID)]=="MS1"){
              MS_filter_qual = FullscanMS1
            } else {
              diff = abs(Peaks$mz[t]-as.numeric(names(Scan_filters)))
              MS_filter_qual = Scan_filters[[match(min(diff),diff)]]
            }
            
            if (T%in%grepl(".mzML", files)) {
              if(MS_filter_qual==FullscanMS1){
                EIC_qual = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                                            grab_what = "EIC",
                                            mz = qual_peaks$mz[q],
                                            ppm = ppm_val,
                                            rtrange = NULL,
                                            prefilter = -1)
                EIC_qual = data.frame("rt"=EIC_qual$EIC$rt,"int"=EIC_qual$EIC$int)
              } else {
                EIC_qual = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                                            grab_what = "MS2",
                                            mz = qual_peaks$mz[q],
                                            ppm = ppm_val,
                                            rtrange = NULL,
                                            prefilter = -1)
                EIC_qual = data.frame("rt"=EIC_qual$MS2$rt,"int"=EIC_qual$MS2$int,"mz"=EIC_qual$MS2$fragmz,"prec"=EIC_qual$MS2$premz)
                precursor = quan_peak$mz
                EIC_qual = EIC_qual[((EIC_qual$prec>=(precursor-precursor_window))&(EIC_qual$prec<=(precursor+precursor_window))),]
                EIC_qual = EIC_qual[(EIC_qual$mz>=RaMS::pmppm(qual_peaks$mz[q],ppm_val)[1]&EIC_qual$mz<=RaMS::pmppm(qual_peaks$mz[q],ppm_val)[2]),]
                EIC_qual = EIC_qual[,c(1:2)]
              }
              
            } else {
              
              Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                               mass = qual_peaks$mz[q],
                                               tol = ppm_val*2,
                                               filter = MS_filter_qual,
                                               type = "xic")
              
              EIC_qual = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
            }
            
            width_qual = abs(qual_peaks$End_RT_level[q]-qual_peaks$Start_RT_level[q])
            
            if(width_qual==0){
              width_qual = abs(sample_peak$End_RT-sample_peak$Start_RT)/2
            }
            
            EIC_qual$int[EIC_qual$int==Inf] = 0
            
            Background_qual_window = qual_peaks$Width[q]*width_factor_background
            if(Background_qual_window==0){
              Background_qual_window=1
            }
            
            Filter_background =  EIC_qual$rt > (Peaks$Start_RT[t]-Background_qual_window-allowed_shift) & EIC_qual$rt < (Peaks$End_RT[t]+Background_qual_window+allowed_shift)
            
            Filter = EIC_qual$rt < (Peaklist_final$Start_RT[quani]-width_qual) | EIC_qual$rt > (Peaklist_final$End_RT[quani]+width_qual)
            
            bad_chroma_qual = F
            
            EIC_qual_background = dplyr::filter(EIC_qual,Filter_background==T)
            EIC_qual = dplyr::filter(EIC_qual,Filter==F)
            if(length(EIC_qual$int)<3|length(EIC_qual_background$int)<3){
              bad_chroma_qual = T
            }
            if(bad_chroma_qual==F){
              plot(EIC_qual$int~EIC_qual$rt,type="p",main=qual_peaks$Compound[q])
            }
            
            if(bad_chroma_qual==F){
              timesplit = ceiling(length(EIC_qual$rt)/width_smoothing)
              EIC_quals = list()
              EIC_quals_smooth = list()
              Start = 1
              for(i in 1:width_smoothing){
                EIC_quals[[i]] = EIC_qual[(Start:(Start+timesplit)),]
                EIC_quals[[i]] = EIC_quals[[i]][!is.na(EIC_quals[[i]]$rt)&!is.na(EIC_quals[[i]]$int),]
                timerange = max(EIC_quals[[i]]$rt)-min(EIC_quals[[i]]$rt)
                bandwidth = 0.1*timerange
                EIC_quals_smooth[[i]] = as.data.frame(ksmooth(EIC_quals[[i]]$rt,EIC_quals[[i]]$int,kernel="normal",bandwidth = bandwidth))
                EIC_quals_smooth[[i]] = EIC_quals_smooth[[i]][!is.na(EIC_quals_smooth[[i]]$x)&!is.na(EIC_quals_smooth[[i]]$y),]
                Start=Start+timesplit
                if(Start>length(EIC_qual$rt)){
                  break
                }
              }
              
              EIC_qual_smooth = dplyr::bind_rows(EIC_quals_smooth)
              colnames(EIC_qual_smooth) = c("rt","int")
              lines(EIC_qual_smooth$rt,EIC_qual_smooth$int,col="pink",type = "p")
              
              threshold = as.numeric(quantile(EIC_qual_background$int,probs=c(normal_background_quantile)))
              for(we in 1:nrow(EIC_qual_background)){
                if(we == 1){
                  if(!is.na(threshold)&EIC_qual_background$int[we]>threshold){
                    sequenci_qual = "b"
                  } else {
                    sequenci_qual = "a"
                  }
                } else {
                  if(!is.na(threshold)&EIC_qual_background$int[we]>threshold){
                    sequenci_qual = append(sequenci_qual,"b")
                  } else {
                    sequenci_qual = append(sequenci_qual,"a")
                  }
                }
              }
              sequencis_qual = paste0(sequenci_qual,collapse = "")
              
              EIC_quals_background = list()
              EIC_quals_background_smooth = list()
              timesplit = ceiling(length(EIC_qual_background$rt)/3)
              Start = 1
              for(i in 1:3){
                EIC_quals_background[[i]] = EIC_qual_background[(Start:(Start+timesplit)),]
                EIC_quals_background[[i]] = EIC_quals_background[[i]][!is.na(EIC_quals_background[[i]]$rt)&!is.na(EIC_quals_background[[i]]$int),]
                timerange = max(EIC_quals_background[[i]]$rt)-min(EIC_quals_background[[i]]$rt)
                bandwidth = 0.1*timerange
                EIC_quals_background_smooth[[i]] = as.data.frame(ksmooth(EIC_quals_background[[i]]$rt,EIC_quals_background[[i]]$int,kernel="normal",bandwidth = bandwidth))
                EIC_quals_background_smooth[[i]] = EIC_quals_background_smooth[[i]][!is.na(EIC_quals_background_smooth[[i]]$x)&!is.na(EIC_quals_background_smooth[[i]]$y),]
                Start=Start+timesplit
                if(Start>length(EIC_qual_background$rt)){
                  break
                }
              }
              
              EIC_qual_background_smooth = dplyr::bind_rows(EIC_quals_background_smooth)
              colnames(EIC_qual_background_smooth) = c("rt","int")
              
              zigzag_positions = as.data.frame(stringr::str_locate_all(sequencis_qual,"aba"))
              zigzag_positions = rbind(zigzag_positions,as.data.frame(stringr::str_locate_all(sequencis_qual,"abba")))
              zigzag_positions = rbind(zigzag_positions,as.data.frame(stringr::str_locate_all(sequencis_qual,"abbba")))
              
              if(nrow(zigzag_positions)>0){
                for(zi in 1:nrow(zigzag_positions)){
                  Start = zigzag_positions$start[zi]
                  End = zigzag_positions$end[zi]
                  sequenci_qual[Start:End] = "c"
                }
              }
              
              abline(v=EIC_qual_background$rt[sequenci_qual=="c"],col=adjustcolor(col="lightgrey",alpha.f = 0.2))
              
              if (T%in%grepl(".mzML", files)){
                mean_Background_qual = as.numeric(quantile(EIC_qual_background$int,na.rm = T,probs = c(normal_background_quantile)))
                max_Background_qual = max(EIC_qual_background$int,na.rm = T)
                if(length(mean_Background_qual)>0&length(max_Background_qual)>0&mean_Background_qual!=0){
                  Background_ratio_qual = max_Background_qual/mean_Background_qual
                } else {
                  Background_ratio_qual = 0
                }
              } else {
                mean_Background_qual = mean(EIC_qual_background$int,na.rm = T)
                max_Background_qual = max(EIC_qual_background$int,na.rm = T)
                
                if(length(mean_Background_qual)>0&length(max_Background_qual)>0&mean_Background_qual!=0){
                  Background_ratio_qual = max_Background_qual/mean_Background_qual
                } else {
                  Background_ratio_qual = 0
                }
              }
              
              if((length(sequenci_qual[sequenci_qual=="c"])>(zigzag_trigger_threshold*length(sequenci_qual)))&(max(EIC_qual_background$int,na.rm = T)!=min(EIC_qual_background$int,na.rm=T))){
                Q_qual = as.numeric(quantile(EIC_qual_background$int[sequenci_qual=="c"],probs=c(higher_background_quantile)))
                Background_qual_signals = EIC_qual_background_smooth$int[EIC_qual_background_smooth$int<Q_qual]
                if(length(Background_qual_signals)==0){
                  Background_qual_signals = EIC_qual_background$int[EIC_qual_background$int<Q_qual]
                }
                Background_qual = sd(Background_qual_signals,na.rm=T)
                Baseline_qual = mean(Background_qual_signals)
              } else if(nrow(EIC_qual_smooth[EIC_qual_smooth$int!=0,])>0){
                if(Background_ratio_qual<3){
                  Q_qual = as.numeric(quantile(EIC_qual_background_smooth$int,probs=c(higher_background_quantile)))
                  if(Q_qual!=0){
                    Background_qual_signals = EIC_qual_background_smooth$int[EIC_qual_background_smooth$int<Q_qual]
                    Background_qual = sd(Background_qual_signals,na.rm=T)
                    Baseline_qual = mean(Background_qual_signals,na.rm=T)
                  } else {
                    Background_qual = 0
                    Baseline_qual = 0
                  }
                  
                } else {
                  Q_qual = as.numeric(quantile(EIC_qual_background_smooth$int,probs=c(0.5)))
                  if(Q_qual!=0){
                    Background_qual_signals = EIC_qual_background_smooth$int[EIC_qual_background_smooth$int<Q_qual]
                    Background_qual = sd(Background_qual_signals,na.rm=T)
                    Baseline_qual = mean(Background_qual_signals,na.rm=T)
                  } else {
                    Background_qual = 0
                    Baseline_qual = 0
                  }
                }
              } else {
                Baseline_qual = 0
                Background_qual = 0
              }
              
              abline(h=Baseline_qual,col="blue",lwd=3)
            }
            if(length(EIC_qual$rt)<3|length(EIC_qual_background$int)<3){
              EIC_qual = data.frame("rt"=c(1,2,3),"int"=c(0,0,0))
              EIC_qual_smooth = EIC_qual
              EIC_qual_background = data.frame("rt"=c(1,2,3),"int"=c(0.000001,0.000001,0.000001))
              EIC_qual_backgroun_smooth = EIC_qual_background
              Background_qual = 0
              Baseline_qual = 0
            }
            
            LOD_qual = as.numeric(Baseline_qual+3*Background_qual)
            LOQ_qual = as.numeric(Baseline_qual+10*Background_qual)
            
            if(bad_chroma_qual==F){
              abline(h=LOD_qual,col="green")
              abline(h=LOQ_qual,col="orange")
            }
            
            guggl = seewave::SAX(EIC_qual$int,alphabet_size = 3,PAA_number = 10)
            
            if(!(T%in%grepl("c",guggl))){
              bad_chroma_qual = F
            }
            
            diffl = abs(EIC_qual_smooth$rt-(min(EIC_qual_smooth$rt,na.rm = T)+0.1))
            density = round(match(min(diffl,na.rm = T),diffl),digits = 0)
            
            #Check if confirming ion is likely and can be expected (depndent on intensity ratios)
            if(LOD_qual!=0){
              if((sample_peak$Height*(qual_peaks$ratio[q]*0.8))<LOQ_qual){
                expected = F
              } else {
                expected = T
              }
            } else {
              if((sample_peak$Height*(qual_peaks$ratio[q]*0.8))<minimum_confirming_peak_height){
                expected = F
              } else {
                expected = T
              }
            }
            
            SAX_result = rep("a",length(EIC_qual$int))
            
            SAX_result_smooth = rep("a",length(EIC_qual_smooth$int))
            
            SAX_result_small = SAX_result
            
            for(l in 1:length(EIC_qual$int)){
              if(EIC_qual$int[l]>LOD_qual){
                SAX_result_small[l] = "b"
              } else {
                SAX_result_small[l] = "a"
              }
            }
            
            if(length(EIC_qual$int)<3|var(EIC_qual$int)==0){
              bad_chroma_qual = T
            }
            
            if(!(T%in%grepl("b",SAX_result_small))){
              bad_chroma_qual = T
            }
            
            if(bad_chroma_qual == F){
              Peaklist = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA)
              
              if(T%in%grepl("b",SAX_result_small)){
                Peak = match(max(EIC_qual_smooth$int,na.rm = T),EIC_qual_smooth$int)
                Peak_rt = EIC_qual_smooth$rt[Peak]
                Peak_int = max(EIC_qual$int[EIC_qual$rt>(Peak_rt-0.05)&EIC_qual$rt<(Peak_rt+0.05)])
                Cutoffs = c()
                if(!is.na(Peak)&Peak_int>LOD_qual){
                  remaining_Peaks = T
                } else {
                  remaining_Peaks = F
                }
                while(remaining_Peaks == T){
                  Peak_start = NA
                  Peak_end = NA
                  Cutoff_r = NA
                  Cutoff_l = NA
                  if(Peak>1){
                    Begin_int = EIC_qual_smooth$int[Peak]
                    Points_int = EIC_qual_smooth$int[(Peak-1):1]
                    Points_index = (Peak-1):1
                    Check = Points_int
                    counter = 0
                    for(c in 1:length(Check)){
                      if(SAX_result_smooth[Points_index[c]]=="b"){
                        if(counter!=0){
                          Peak_start = Points_index[c-counter]
                          break
                        }
                        Peak_start = Points_index[[c]]+1
                        break
                      }
                      if(Check[c]>=(Begin_int*0.7)){
                        SAX_result_smooth[Points_index[c]] = "b"
                        while(counter!=0){
                          SAX_result_smooth[Points_index[c-counter]] = "b"
                          counter = counter-1
                        }
                      } else {
                        counter = counter+1
                      }
                      if(counter==7){
                        Peak_start = Points_index[c-7]
                        if(c==7){
                          Peak_start = Peak
                        }
                        break
                      }
                      if(c==length(Points_index)){
                        if(counter==0){
                          Peak_start = Points_index[length(Points_index)]
                        } else {
                          if(c==counter){
                            Peak_start = 1
                          } else {
                            Peak_start = Points_index[length(Points_index)-counter]
                          }
                        }
                      }
                    }
                  } else {
                    Peak_start = 1
                  }
                  
                  if(Peak < length(SAX_result_smooth)){
                    End_int = EIC_qual_smooth$int[Peak]
                    Points_int = EIC_qual_smooth$int[(Peak+1):length(SAX_result_smooth)]
                    Points_index = (Peak+1):length(SAX_result_smooth)
                    Check = Points_int
                    counter = 0
                    for(c in 1:length(Check)){
                      if(SAX_result_smooth[Points_index[c]]=="b"){
                        if(counter!=0){
                          Peak_end = Points_index[c-counter]
                          break
                        }
                        Peak_end = Points_index[[c]]-1
                        break
                      }
                      if(Check[c]>=(0.7*End_int)){
                        SAX_result_smooth[Points_index[c]] = "b"
                        SAX_result_smooth[Points_index[c-1]] = "b"
                        while(counter!=0){
                          SAX_result_smooth[Points_index[c-counter]] = "b"
                          counter = counter-1
                        }
                      } else {
                        counter = counter+1
                      }
                      if(counter==7){
                        Peak_end = Points_index[c-7]
                        if(c==7){
                          Peak_end = Peak
                        }
                        break
                      }
                      if(c==length(Points_index)){
                        if(counter==0){
                          Peak_end = Points_index[length(Points_index)]
                        } else {
                          if(c==counter){
                            Peak_end = length(SAX_result_smooth)
                          } else {
                            Peak_end = Points_index[length(Points_index)-counter]
                          }
                        }
                      }
                    }
                  } else {
                    Peak_end = length(SAX_result_smooth)
                  }
                  
                  SAX_result_smooth[Peak] = "b"
                  
                  if(Peak_start > 1){
                    Begin_int = EIC_qual_smooth$int[Peak_start]
                    Points_int = EIC_qual_smooth$int[(Peak_start):1]
                    Points_index = (Peak_start):1
                    Check = Points_int
                    counter = 0 #counts negative trends
                    counter2 = 0 #counts positive trends
                    threshold = NA
                    for(c in 2:length(Check)){
                      if(SAX_result_smooth[Points_index[c]]=="b"){
                        if(counter!=0){
                          Cutoff_l = Points_index[c-counter]
                          break
                        }
                        Cutoff_l = Points_index[c-1]
                        break
                      }
                      if(counter==0){
                        if(Check[c]<(Check[c-1])){
                          SAX_result_smooth[Points_index[c]] = "b"
                          SAX_result_smooth[Points_index[c-1]] = "b"
                          Cutoff_l = Points_index[c]
                        } else {
                          counter = counter+1
                          threshold = Points_int[c-1]*0.6
                          Cutoff_l = Points_index[c-1]
                        }
                      }else{
                        if(Check[c]>(1.4*1.6666666*threshold)){
                          counter2 = counter2+1
                        }
                        if(counter2>(0.25*density+1)){
                          threshold=0
                        }
                        if(Check[c]<threshold){
                          SAX_result_smooth[Points_index[c]] = "b"
                          SAX_result_smooth[Points_index[c-1]] = "b"
                          Cutoff_l = Points_index[c]
                          threshold = NA
                          counter2 = 0
                          while(counter!=0){
                            SAX_result_smooth[Points_index[c-counter]] = "b"
                            counter = counter-1
                          }
                        } else {
                          counter = counter+1
                        }
                        if(counter==(5*density+1)){
                          Cutoff_l = Points_index[(c-((5*density+1)-1))]
                          break
                        }
                      }
                    }
                    if(is.na(Cutoff_l)){
                      Cutoff_l = 1
                    }
                  } else {
                    Cutoff_l = 1
                  }
                  
                  if(Peak_end < length(SAX_result_smooth)){
                    End_int = EIC_qual_smooth$int[Peak_end]
                    Points_int = EIC_qual_smooth$int[Peak_end:length(SAX_result_smooth)]
                    Points_index = Peak_end:length(SAX_result_smooth)
                    Check = Points_int
                    counter = 0 #counts negative trends
                    counter2 = 0 #counts positive trends
                    threshold = NA
                    for(c in 2:length(Check)){
                      if(SAX_result_smooth[Points_index[c]]=="b"){
                        if(counter!=0){
                          Cutoff_r = Points_index[c-counter]
                          break
                        }
                        Cutoff_r = Points_index[c-1]
                        break
                      }
                      if(counter==0){
                        if(Check[c]<(Check[c-1])){
                          SAX_result_smooth[Points_index[c]] = "b"
                          SAX_result_smooth[Points_index[c-1]] = "b"
                          Cutoff_r = Points_index[c]
                        } else {
                          counter = counter+1
                          if(is.na(threshold)){
                            threshold = Points_int[c-1]*0.6
                          }
                          Cutoff_r = Points_index[c-1]
                        }
                      }else{
                        if(Check[c]>(1.4*1.6666666*threshold)){
                          counter2 = counter2+1
                        }
                        if(counter2>(0.25*density+1)){
                          threshold=0
                        }
                        if(Check[c]<threshold){
                          SAX_result_smooth[Points_index[c]] = "b" 
                          SAX_result_smooth[Points_index[c-1]] = "b"
                          Cutoff_r = Points_index[c]
                          threshold = NA
                          counter2 = 0
                          while(counter!=0){
                            SAX_result_smooth[Points_index[c-counter]] = "b"
                            counter = counter-1
                          }
                        } else {
                          counter = counter+1
                        }
                        if(counter==(5*density+1)){
                          Cutoff_r = Points_index[(c-((5*density+1)-1))]
                          break
                        }
                      }
                    }
                    if(is.na(Cutoff_r)){
                      Cutoff_r = length(EIC_qual_smooth$rt)
                    }
                  } else {
                    Cutoff_r = length(EIC_qual_smooth$rt)
                  }
                  
                  Pairs = length(Cutoffs)/2
                  if(Pairs>0){
                    for(pa in 1:Pairs){
                      pairis = data.frame("Start"=seq(1,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2),"End"=seq(2,plyr::round_any(Pairs*2,accuracy = 100,f=ceiling),by=2))
                      Range = Cutoffs[pairis$Start[pa]]:Cutoffs[pairis$End[pa]]
                      if(Cutoff_l%in%Range){
                        Cutoff_l = max(Range)
                      }
                      if(Cutoff_r%in%Range){
                        Cutoff_r = min(Range)
                      }
                    }
                  }
                  
                  Cutoffs = append(Cutoffs,c(Cutoff_l,Cutoff_r))
                  
                  Indizes = 1:length(EIC_qual_smooth$int)
                  
                  Filter = SAX_result_smooth!="b"
                  
                  Indizes = Indizes[Filter]
                  
                  Peak = Indizes[match(max(EIC_qual_smooth$int[SAX_result_smooth!="b"]),EIC_qual_smooth$int[Filter])]
                  
                  Peak_rt = EIC_qual_smooth$rt[Peak]
                  
                  Peak_int = max(EIC_qual$int[EIC_qual$rt>(Peak_rt-0.05)&EIC_qual$rt<(Peak_rt+0.05)])
                  
                  if(!is.na(Peak)&Peak_int>LOD_qual){
                    remaining_Peaks = T
                  } else {
                    remaining_Peaks = F
                  }
                }
                SAX_result_smooth[Cutoffs] = "a"
                if(1%in%Cutoffs){
                  SAX_result_smooth[1] = "b"
                }
                if(length(SAX_result_smooth)%in%Cutoffs){
                  SAX_result_smooth[length(SAX_result_smooth)] = "b"
                }
              }
              
              Start = 1
              
              while(Start<length(SAX_result_smooth)){
                peakbegin = match("b",SAX_result_smooth[Start:length(SAX_result_smooth)])+Start-1
                if(is.na(peakbegin)){
                  break
                }
                peakend = match("a",SAX_result_smooth[peakbegin:length(SAX_result_smooth)])+peakbegin-2
                if(is.na(peakend)){
                  peakend = length(SAX_result_smooth)
                  if(peakbegin==peakend){
                    break
                  }
                }
                if(peakbegin==peakend){
                  Start = peakend+1
                  next
                }
                Starttime = EIC_qual_smooth$rt[peakbegin]
                Endtime = EIC_qual_smooth$rt[peakend]
                SAX_result[EIC_qual$rt>=Starttime&EIC_qual$rt<=Endtime] = "b"
                Start = peakend + 1
              }
              
              SAX_result[EIC_qual$int<=extended_baseline_factor*Baseline_qual] = "a"
              
              Start = 1
              while(Start<length(SAX_result)){
                peakbegin = match("b",SAX_result[Start:length(SAX_result)])+Start-1
                if(is.na(peakbegin)){
                  break
                }
                peakend = match("a",SAX_result[peakbegin:length(SAX_result)])+peakbegin-2
                if(is.na(peakend)){
                  peakend = length(SAX_result)
                  if(peakbegin==peakend){
                    break
                  }
                }
                if(peakbegin==peakend){
                  Start = peakend+1
                  next
                }
                if(peakbegin>1){
                  if(EIC_qual$int[peakbegin-1]<=(extended_baseline_factor*Baseline_qual)){
                    peakbegin = peakbegin - 1
                  }
                }
                if(peakend<length(SAX_result)){
                  if(EIC_qual$int[peakend+1]<=(extended_baseline_factor*Baseline_qual)){
                    peakend = peakend + 1
                  }
                }
                Starttime = EIC_qual$rt[peakbegin]
                Endtime = EIC_qual$rt[peakend]
                SAX_result[EIC_qual$rt>=Starttime&EIC_qual$rt<=Endtime] = "b"
                Start = peakend + 1
              }
              
              Cutoffs_time = EIC_qual_smooth$rt[Cutoffs]
              
              for(klo in 1:length(Cutoffs_time)){
                diff = abs(EIC_qual$rt-Cutoffs_time[klo])
                time = EIC_qual$rt[match(min(diff,na.rm = T),diff)]
                SAX_result[EIC_qual$rt==time] = "a"
              }
              
              #Test whether it is plausible that there are traces of analyte or just non-specific nearby peaks
              if(is.null(Cutoffs)){
                Cutoffs = NA
              }
              if(!is.na(Cutoffs)&length(Cutoffs)>0){
                n_maxima = length(Cutoffs)/2
                cutoff_maxima = rep(NA,n_maxima)
                cutoff_indices = matrix(c(1:length(Cutoffs)),ncol=2,byrow = T)
                for(c in 1:n_maxima){
                  cutoff_maxima[c] = max(EIC_qual_smooth$int[Cutoffs[cutoff_indices[c,1]]:Cutoffs[cutoff_indices[c,2]]])
                }
                cutoff_maxima_rt = EIC_qual_smooth$rt[which(EIC_qual_smooth$int%in%cutoff_maxima)]
                integrate_peak_rt = F
                for(d in 1:n_maxima){
                  if(allowed_shift>quan_peak$RT_tol){
                    cutoff_rt = ifelse(cutoff_maxima_rt[d]<=quan_peak$RT+allowed_shift&cutoff_maxima_rt[d]>=quan_peak$RT-allowed_shift,cutoff_maxima_rt[d],NA)
                  } else {
                    cutoff_rt = ifelse(cutoff_maxima_rt[d]<=quan_peak$RT+quan_peak$RT_tol&cutoff_maxima_rt[d]>=quan_peak$RT-quan_peak$RT_tol,cutoff_maxima_rt[d],NA)
                  }
                  if(is.na(cutoff_rt)){
                    next
                  }
                  diff_cutoff = abs(cutoff_rt-cutoff_maxima_rt)
                  diff_cutoff = abs(cutoff_rt-cutoff_maxima_rt)
                  if((min(diff_cutoff[diff_cutoff!=0])<quan_peak$Width*width_factor_background)){
                    integrate_peak_rt = T
                    break
                  } 
                }
              } else {
                integrate_peak_rt = F
              }
              if(integrate_peak_rt){
                if(quan_peak$Start_RT_level==quan_peak$End_RT_level){
                  EIC_qual_sub = EIC_qual[EIC_qual$rt>=quan_peak$Start_RT&EIC_qual$rt<=quan_peak$End_RT,]
                  EIC_qual_sub_narrow = EIC_qual[EIC_qual$rt>=quan_peak$Start_RT_level&EIC_qual$rt<=quan_peak$End_RT_level,]
                  EIC_qual_not_peak = EIC_qual_smooth[EIC_qual_smooth$rt%in%cutoff_maxima_rt[!(cutoff_maxima_rt>min(EIC_qual_sub$rt)&cutoff_maxima_rt<max(EIC_qual_sub$rt))],]
                  if(nrow(EIC_qual_not_peak)>0){
                    cutoff_background_ratio = quantile(EIC_qual_sub$int,probs = c(0.95),na.rm=T)/quantile(EIC_qual_not_peak$int,probs = c(0.95),na.rm=T)
                    if(is.na(cutoff_background_ratio)|is.nan(cutoff_background_ratio)){
                      cutoff_background_ratio = 0
                    }
                  } else {
                    cutoff_background_ratio = minimum_background_ratio
                  }
                  intensities_in_centre = nrow(EIC_qual_sub_narrow[EIC_qual_sub_narrow$int>LOD,])>0
                  if((cutoff_background_ratio<minimum_background_ratio)&intensities_in_centre){
                    SAX_aimed_peak = rep("a",nrow(EIC_qual_sub))
                  } else {
                    SAX_aimed_peak = SAX_result[EIC_qual$rt>=quan_peak$Start_RT&EIC_qual$rt<=quan_peak$End_RT]
                  }
                } else {
                  EIC_qual_sub = EIC_qual[EIC_qual$rt>=quan_peak$Start_RT&EIC_qual$rt<=quan_peak$End_RT,]
                  EIC_qual_sub_narrow = EIC_qual[EIC_qual$rt>=quan_peak$Start_RT_level&EIC_qual$rt<=quan_peak$End_RT_level,]
                  EIC_qual_not_peak = EIC_qual_smooth[EIC_qual_smooth$rt%in%cutoff_maxima_rt[!(cutoff_maxima_rt>min(EIC_qual_sub$rt)&cutoff_maxima_rt<max(EIC_qual_sub$rt))],]
                  if(nrow(EIC_qual_not_peak)>0){
                    cutoff_background_ratio = quantile(EIC_qual_sub$int,probs = c(0.95),na.rm=T)/quantile(EIC_qual_not_peak$int,probs = c(0.95),na.rm=T)
                    if(is.na(cutoff_background_ratio)|is.nan(cutoff_background_ratio)){
                      cutoff_background_ratio = 0
                    }
                  } else {
                    cutoff_background_ratio = minimum_background_ratio
                  }
                  intensities_in_centre = nrow(EIC_qual_sub_narrow[EIC_qual_sub_narrow$int>LOD,])>0
                  if((cutoff_background_ratio<minimum_background_ratio)&intensities_in_centre){
                    SAX_aimed_peak = rep("a",nrow(EIC_qual_sub))
                  } else {
                    SAX_aimed_peak = SAX_result[EIC_qual$rt>=quan_peak$Start_RT_level&EIC_qual$rt<=quan_peak$End_RT_level]
                  }
                }
                
                if(T%in%grepl("b",SAX_aimed_peak)){
                  if(quan_peak$Start_RT_level==quan_peak$End_RT_level){
                    Start = match("b",SAX_aimed_peak)
                    End = length(SAX_aimed_peak)-match("b",rev(SAX_aimed_peak))+1
                    SAX_aimed_peak[Start:End] = "b"
                  } else {
                    SAX_aimed_peak[1:length(SAX_aimed_peak)] = "b"
                  }
                  SAX_result[EIC_qual$rt>=quan_peak$Start_RT&EIC_qual$rt<=quan_peak$End_RT] = SAX_aimed_peak
                }
              }
              
              #If two peaks are separated and Cutoff > LOQ, apply it 
              Cutoffs_time = EIC_qual_smooth$rt[Cutoffs]
              Cutoffs_intensity = EIC_qual_smooth$int[Cutoffs]
              if(LOQ_qual>0){
                Cutoffs_time = Cutoffs_time[Cutoffs_intensity>LOQ_qual]
              } else {
                Cutoffs_time = Cutoffs_time[Cutoffs_intensity>minimum_cutoff_intensity_factor*quan_peak$Height]
              }
              
              
              for(klo in 1:length(Cutoffs_time)){
                diff = abs(EIC_qual$rt-Cutoffs_time[klo])
                time = EIC_qual$rt[match(min(diff,na.rm = T),diff)]
                SAX_result[EIC_qual$rt==time] = "a"
              }
              
              Start = 1
              End = length(SAX_result)
              Count = 1
              Already = 0
              
              while(5==5){
                not_a = match(letters[2:26],SAX_result[Start:End])+Already
                k = not_a[order(not_a,decreasing = F)[1]]
                k = k
                if(!is.na(k)){
                  k = k
                }
                l = match("a",SAX_result[k+1:End])+k-Start+Already+1
                if(!is.na(l)){
                  l = l
                }
                #test for gaps
                m = k - 1
                if(is.na(m)){
                  break
                }
                if(m==0){
                  m = 1
                }
                if(abs(EIC_qual$rt[k]-EIC_qual$rt[m])>0.01){
                  m = k
                }
                if(is.na(l)){
                  if(mean(!is.na(not_a))!=0){
                    level = quan_peak$level
                    l = End
                    Peaklist$Start_RT[Count] = EIC_qual$rt[m]
                    Peaklist$End_RT[Count] = EIC_qual$rt[l]
                    Smoothie = EIC_qual_smooth[EIC_qual_smooth$rt>=EIC_qual$rt[m]&EIC_qual_smooth$rt<=EIC_qual$rt[l],]
                    Peakl = EIC_qual[EIC_qual$rt>=EIC_qual$rt[m]&EIC_qual$rt<=EIC_qual$rt[l],]
                    Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
                    Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
                    if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
                      Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
                    } else {
                      Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
                    }
                    Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                                        error=function(e){"too_small"})
                    Peaklist$Nr_of_Points[Count] = length(EIC_qual$int[m:l])
                    Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
                    Peaklist$Height[Count] = max(EIC_qual$int[m:l])
                    Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
                    Peaklist$RT[Count] = median(Peaktop$rt)
                    Smoothie1 = Smoothie
                    Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
                    #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
                    tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                                             resp = Smoothie1$int,
                                             bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
                    
                    RMSD_Ratio = Background_qual/Baseline_qual*100
                    if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
                      RMSD_Ratio = 0
                    }
                    if(RMSD_Ratio<10){
                      Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
                    } else {
                      Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
                    }
                    
                    if(is.na(Ratio)){
                      Ratio = 0
                    }
                    Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
                    Threshold_Differenci = plyr::round_any(Differenci,0.01,f=floor)
                    Threshold_Ratio = plyr::round_any(Ratio,0.1,f=ceiling)
                    Threshold_Accepted = Threshold_Differenci<=0.1&Threshold_Ratio>=2
                    if(grepl("GC",sample)){
                      if((Differenci<=0.01&Ratio<3)|Threshold_Accepted){
                        Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
                      }
                    } else {
                      if((Differenci<=0.05&Ratio<3)|Threshold_Accepted){
                        Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
                      }
                    }
                    
                    #Constant = no real peak, too flat
                    Constant_AIC = tcpl_fit$cnst_aic
                    #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
                    Gainloss_AIC = tcpl_fit$gnls_aic
                    
                    if(!is.na(Gainloss_AIC)){
                      if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
                        Gainloss_AIC = NA
                      }
                    }
                    id = m:l
                    Peaklist$Area[Count] = sum(diff(EIC_qual$rt[id])*rollmean(EIC_qual$int[id],2))
                    Count = Count+1
                    Peaklist[Count,] = NA
                    if(Ratio<3|is.na(Gainloss_AIC)){
                      Peaklist[(Count-1),] = NA
                    }
                    break
                  } else{
                    break
                  }
                }
                level = quan_peak$level
                Peaklist$Start_RT[Count] = EIC_qual$rt[m]
                Peaklist$End_RT[Count] = EIC_qual$rt[l]
                Smoothie = EIC_qual_smooth[EIC_qual_smooth$rt>=EIC_qual$rt[m]&EIC_qual_smooth$rt<=EIC_qual$rt[l],]
                Peakl = EIC_qual[EIC_qual$rt>=EIC_qual$rt[m]&EIC_qual$rt<=EIC_qual$rt[l],]
                Peaklist$Start_RT_level[Count] = Peakl$rt[match(T,Peakl$int>=level*max(Peakl$int,na.rm = T))]
                Peaklist$End_RT_level[Count] = Peakl$rt[length(Peakl$int)-match(T,Peakl$int[length(Peakl$int):1]>=level*max(Peakl$int,na.rm = T))+1]
                if(Peaklist$Start_RT_level[Count]==Peaklist$End_RT_level[Count]){
                  Peaki = EIC[EIC$rt>=Peaklist$Start_RT[Count]&EIC$rt<=Peaklist$End_RT[Count],]
                } else {
                  Peaki = EIC[EIC$rt>=Peaklist$Start_RT_level[Count]&EIC$rt<=Peaklist$End_RT_level[Count],]
                }
                Peaklist$Sequence[Count] = tryCatch(paste0(seewave::SAX(relative(Peaki$int),alphabet_size=alphabet_size,PAA_number = (alphabet_size+1),breakpoints="quantiles"),collapse = ""),
                                                    error=function(e){"too_small"})
                Peaktop = Peakl[Peakl$int>=0.8*max(Peakl$int,na.rm = T),]
                Peaklist$RT[Count] = median(Peaktop$rt)
                
                Smoothie1 = Smoothie
                Smoothie1$int = ((Smoothie1$int - min(Smoothie1$int,na.rm=T)) / (max(Smoothie1$int,na.rm=T) - min(Smoothie1$int,na.rm=T))) * 100
                #tcpl needs 0 - 100 % and a log-logistic range; ToxCast used 0 - 100 ???M
                tcpl_fit = tcpl::tcplFit(logc = seq(-1, 2, length.out = nrow(Smoothie1)),
                                         resp = Smoothie1$int,
                                         bmad = quantile(Smoothie1$int,probs = c(0.1),na.rm = T))
                
                RMSD_Ratio = Background_qual/Baseline_qual*100
                if(is.na(RMSD_Ratio)|is.nan(RMSD_Ratio)){
                  RMSD_Ratio = 0
                }
                if(RMSD_Ratio<10){
                  Ratio = ceiling(max(Smoothie$int,na.rm=T)/as.numeric(quantile(Smoothie$int,probs = c(0.01))))
                } else {
                  Ratio = (as.numeric(quantile(Smoothie$int,probs = c(0.8)))/as.numeric(quantile(Smoothie$int,probs = c(0.1))))
                }
                
                if(is.na(Ratio)){
                  Ratio = 0
                }
                Differenci = abs(Peaks$RT[t]-Peaklist$RT[Count])
                Threshold_Differenci = plyr::round_any(Differenci,0.01,f=floor)
                Threshold_Ratio = plyr::round_any(Ratio,0.1,f=ceiling)
                Threshold_Accepted = Threshold_Differenci<=0.1&Threshold_Ratio>=1.5
                if(grepl("GC",sample)){
                  if((Differenci<=0.01&Ratio<3)|Threshold_Accepted){
                    Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
                  }
                } else {
                  if((Differenci<=0.05&Ratio<3)|Threshold_Accepted){
                    Ratio = 3 #if you have a Peak exactly when you would expect it, keep this peak so that neighboring peaks are less likely to be selected
                  }
                }
                
                #Constant = no real peak, too flat
                Constant_AIC = tcpl_fit$cnst_aic
                #Gainloss = Increasing, then decreasing (in significant manner) -> likely peak
                Gainloss_AIC = tcpl_fit$gnls_aic
                if(!is.na(Gainloss_AIC)){
                  if(abs(Gainloss_AIC-Constant_AIC)<5|Gainloss_AIC>Constant_AIC){
                    Gainloss_AIC = NA
                  }
                }
                #If gnls_aic = NA or is bigger than Constant_AIC, there is likely NO valid peak
                Peaklist$Nr_of_Points[Count] = length(EIC_qual$int[m:l])
                Peaklist$Width[Count] = Peaklist$End_RT_level[Count]-Peaklist$Start_RT_level[Count]
                Peaklist$Height[Count] = max(EIC_qual$int[m:l])
                id = m:l
                Peaklist$Area[Count] = sum(diff(EIC_qual$rt[id])*rollmean(EIC_qual$int[id],2))
                Start = l+1
                Count = Count+1
                Already = l
                Peaklist[Count,] = NA
                if(Ratio<3|is.na(Gainloss_AIC)){
                  Peaklist[(Count-1),] = NA
                }
              }
              Peaklist = Peaklist[!is.na(Peaklist$RT)&Peaklist$Sequence!="too_small",]
              
              if(ID%in%Multiple$ID&nrow(Peaklist)>0){
                if(Multiple$multiple[match(ID,Multiple$ID)]==0){
                  for(m in nrow(Peaklist)){
                    Peaklist$RT[m] = median(c(Peaklist$Start_RT_level[m],Peaklist$End_RT_level[m]),na.rm = T)
                  }
                }
              }
              
              if(nrow(Peaklist)==0){
                Peaklist[1,] = NA
              }
              
              Check_sequence = (stringr::str_count(Peaklist$Sequence,"a")<maximum_nr_of_a|sum_letters(Peaklist$Sequence,alphabet_size)>=minimum_nr_of_high_intensity_letters)&Peaklist$Area > minimum_peak_area
              
              Filter = !is.na(Peaklist$Area) & Peaklist$Width < maximum_peak_width & Check_sequence
              Peaklist_final_qual = filter(Peaklist,Filter==T)
              if(mean(Filter)==0 | length(Filter)==0){
                Peaklist_final_qual = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
              }
            } else {
              Peaklist_final_qual = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
            }
            Peaklist_final_qual$Compound = qual_peaks$Compound[q]
            Peaklist_final_qual$mz = qual_peaks$mz[q]
            if(!is.na(Peaklist_final_qual$Height[1])){
              if(mean(Filter)==0 | length(Filter)==0){
                Peaklist_final_qual = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA)
                Peaklist_final_qual$Compound = qual_peaks$Compound[q]
                Peaklist_final_qual$mz = qual_peaks$mz[q]
              }
            }
            
            Peaklist_final_qual = Peaklist_final_qual[!is.na(Peaklist_final_qual$RT),]
            
            if(length(Peaklist_final_qual$Comment)>1){
              ID = as.numeric(substr(Peaks$Compound[t],1,4))
              if(is.na(ID)){ #Internal Standards
                ID = substr(Peaks$Compound[t],1,4)
              }
              if(ID%in%Multiple$ID){
                if(Multiple$multiple[match(ID,Multiple$ID)]==0){
                  #Since the loop is for each peak, no merging if multiple ones are relevant
                  
                  #first check if very small peak is included by mistake
                  #Peaklist_final_qual = Peaklist_final_qual[Peaklist_final_qual$Height>inner_intensity_ratio_multiple_peaks*max(Peaklist_final_qual$Height,na.rm=T),]
                  #potential_width = abs(Peaklist_final_qual$End_RT_level[nrow(Peaklist_final_qual)]-Peaklist_final_qual$Start_RT_level[1])
                  #aimed_width = quan_peak$End_RT_level-quan_peak$Start_RT_level
                  
                  #while(potential_width/aimed_width>1.2&nrow(Peaklist_final_qual)>1){
                  #  #Remove the peak which is most far away from center if widths are deviating (likely unspecific peak present)
                  #  shift_corrected_times = Peaks$RT[t]+Peaklist_final_qual$intensity_shift+IS_dependent_shift_compound
                  #  center = shift_corrected_times[1]
                  #  diffs = abs(Peaklist_final_qual$RT-center)
                  #  Peaklist_final_qual = Peaklist_final_qual[diffs!=max(diffs,na.rm = T),]
                  #  potential_width = abs(Peaklist_final_qual$End_RT_level[nrow(Peaklist_final_qual)]-Peaklist_final_qual$Start_RT_level[1])
                  #}
                  diff = abs(Peaklist_final_qual$RT-sample_peak$RT)
                  Peaklist_final_qual = Peaklist_final_qual[match(min(diff),diff),]
                  quali_peak = Peaklist_final_qual
                  #Since the loop is for each peak, no merging if multiple ones are relevant
                  #quali_peak = Peaklist_final_qual[1,]
                  #quali_peak$Start_RT = min(Peaklist_final_qual$Start_RT,na.rm = T)
                  #quali_peak$RT = median(Peaklist_final_qual$RT,na.rm = T)
                  #quali_peak$End_RT = max(Peaklist_final_qual$End_RT,na.rm = T)
                  #quali_peak$Sequence = paste0(Peaklist_final_qual$Sequence,collapse = "")
                  #quali_peak$Width = quali_peak$End_RT-quali_peak$Start_RT
                  #quali_peak$Height = max(Peaklist_final_qual$Height,na.rm = T)
                  #quali_peak$Area = sum(Peaklist_final_qual$Area,na.rm = T)
                } else {
                  diff = abs(Peaklist_final_qual$RT-Peaks$RT[t])
                  Peaklist_final_qual = Peaklist_final_qual[match(min(diff),diff),]
                  quali_peak = Peaklist_final_qual
                }
              } else {
                diff = abs(Peaklist_final_qual$RT-Peaks$RT[t])
                Peaklist_final_qual = Peaklist_final_qual[match(min(diff),diff),]
                quali_peak = Peaklist_final_qual
              }
            } else {
              quali_peak = Peaklist_final_qual
            }
            
            quali_peak = quali_peak[!is.na(quali_peak$RT),]
            
            if(length(quali_peak$RT)>0){
              diff = abs(quali_peak$RT-sample_peak$RT)
            } else {
              diff = Inf
            }
            
            if(bad_chroma_qual==F&nrow(quali_peak[!is.na(quali_peak$Height),])>0){
              for(qlk in 1:nrow(quali_peak)){
                abline(v=quali_peak$Start_RT[qlk],col="purple")
                abline(v=quali_peak$End_RT[qlk],col="purple")
                Peak = EIC_qual[EIC_qual$rt>=quali_peak$Start_RT[qlk]&EIC_qual$rt<=quali_peak$End_RT[qlk],]
                lines(Peak$rt,Peak$int,type="l",col="purple",lwd=2)
              }
            }
            
            if(length(quali_peak$RT)==0&expected==T){
              Peaklist_final$confirmed[quani] = -Inf
              break
            }
            
            if(ID%in%Multiple$ID){
              k = match(ID,Multiple$ID)
              if(Multiple$multiple[k]==0){
                allowed_diff =  0.1
              } else {
                if(T%in%grepl("GC",sample)){
                  allowed_diff = 0.025
                } else {
                  allowed_diff = 0.1
                }
              }
            } else {
              if(T%in%grepl("GC",sample)){
                allowed_diff = 0.025
              } else {
                allowed_diff = 0.1
              }
            }
            
            if(diff>allowed_diff&expected==T){
              Peaklist_final$confirmed[quani] = -Inf
              break
            }
            
            if((diff<allowed_diff&expected==T)|(expected==F&diff<allowed_diff)){
              Peaklist_final$confirmed[quani] = Peaklist_final$confirmed[quani]+1
            }
          }
        }
      }
    }
  }
  #If there are confirming ions present, this is one major argument for peak x
  Peaklist_final = Peaklist_final[!is.na(Peaklist_final$confirmed),]
  if(nrow(Peaklist_final)==0){
    Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
    Peaklist_final$Compound = quan_peak$Compound
    Peaklist_final$mz = quan_peak$mz
    Peaklist_final$Comment = paste0(Peaklist_final$Comment," expected confirming ion(s) not found;")
  }
  if(!is.na(Peaklist_final$RT)){
    RT_diffis = abs(Peaklist_final$RT-Peaks$RT[t])
    if(T%in%grepl("GC",sample)){
      perfect_RT_match = T%in%RT_diffis<=0.05
    } else {
      perfect_RT_match = T%in%RT_diffis<=0.1
    }
    perfect_match_confirmed = match(min(RT_diffis,na.rm=T),RT_diffis)==match(max(Peaklist_final$confirmed),Peaklist_final$confirmed)
    if(!perfect_RT_match|perfect_match_confirmed){
      Peaklist_final = Peaklist_final[Peaklist_final$confirmed==max(Peaklist_final$confirmed,na.rm=T),]
      if(nrow(Peaklist_final)==0){
        Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
        Peaklist_final$Compound = quan_peak$Compound
        Peaklist_final$mz = quan_peak$mz
        Peaklist_final$Comment = paste0(Peaklist_final$Comment," expected confirming ion(s) not found;")
      }
      Peaklist_final = Peaklist_final[!is.na(Peaklist_final$RT),]
      if(nrow(Peaklist_final)==0){
        Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
        Peaklist_final$Compound = quan_peak$Compound
        Peaklist_final$mz = quan_peak$mz
      }
    }
  }
  
  #Create final sample peak
  if(length(Peaklist_final$Comment)>1){
    ID = as.numeric(substr(Peaks$Compound[t],1,4))
    if(is.na(ID)){
      ID = paste0("IS",substr(Peaks$Compound[t],3,4))
    }
    if(ID%in%Multiple$ID){
      if(Multiple$multiple[match(ID,Multiple$ID)]==0){
        shift_corrected_times = Peaks$RT[t]+Peaklist_final$intensity_shift+IS_dependent_shift_compound
        Peaklist_final = Peaklist_final[Peaklist_final$RT>(shift_corrected_times-RT_range)&Peaklist_final$RT<(shift_corrected_times+RT_range),]
        if(nrow(Peaklist_final)==0){
          Peaklist_final = data.frame("Compound"=NA,"mz"=NA, "Comment" = NA, "Start_RT" = NA,"RT" = NA,"End_RT" = NA,"Start_RT_level"=NA,"End_RT_level"=NA,"Sequence" = NA,"Nr_of_Points"=NA,"Width" = NA,"Height" = NA,"Area" = NA,"intensity_shift"=NA, "LOD" = NA, "LOQ" = NA,"confirmed"=NA)
          Peaklist_final$Compound = quan_peak$Compound
          Peaklist_final$mz = quan_peak$mz
        } else {
          #first check if very small peak is included by mistake
          
          Peaklist_final = Peaklist_final[Peaklist_final$Height>inner_intensity_ratio_multiple_peaks*max(Peaklist_final$Height,na.rm=T),]
      
          potential_width = abs(Peaklist_final$End_RT_level[nrow(Peaklist_final)]-Peaklist_final$Start_RT_level[1])
          aimed_width = quan_peak$End_RT_level-quan_peak$Start_RT_level
          
          while(potential_width/aimed_width>1.2&nrow(Peaklist_final)>1){
            #Remove the peak which is most far away from center if widths are deviating (likely unspecific peak present)
            shift_corrected_times = Peaks$RT[t]+Peaklist_final$intensity_shift+IS_dependent_shift_compound
            center = shift_corrected_times[1]
            diffs = abs(Peaklist_final$RT-center)
            Peaklist_final = Peaklist_final[diffs!=max(diffs,na.rm = T),]
            potential_width = abs(Peaklist_final$End_RT_level[nrow(Peaklist_final)]-Peaklist_final$Start_RT_level[1])
          }
          sample_peak = Peaklist_final[1,]
          sample_peak$Start_RT = min(Peaklist_final$Start_RT,na.rm = T)
          sample_peak$Start_RT_level = min(Peaklist_final$Start_RT_level,na.rm=T)
          sample_peak$RT = median(Peaklist_final$RT,na.rm = T)
          sample_peak$End_RT = max(Peaklist_final$End_RT,na.rm = T)
          sample_peak$End_RT_level = max(Peaklist_final$End_RT_level,na.rm=T)
          sample_peak$Sequence = SAX_consensus(Peaklist_final$Sequence)
          sample_peak$Comment = Peaklist_final$Comment[1]
          sample_peak$Width = sample_peak$End_RT-sample_peak$Start_RT
          sample_peak$Height = max(Peaklist_final$Height,na.rm = T)
          sample_peak$Area = sum(Peaklist_final$Area,na.rm = T)
        }
      } else {
        shift_corrected_times = Peaks$RT[t]+Peaklist_final$intensity_shift+IS_dependent_shift_compound
        diff = abs(Peaklist_final$RT-shift_corrected_times)
        diff.grid = expand.grid(diff,diff)
        diff.grid$diff = abs(diff.grid$Var1-diff.grid$Var2)
        diff.grid = diff.grid[diff.grid$diff!=0,]
        differences = unique(diff.grid$diff)
        nearest_peaks = min(diff)==diff.grid$Var1[diff.grid$diff==min(diff.grid$diff)]
        if(min(differences)<0.05&T%in%nearest_peaks){
          minidist = rep(NA,nrow(Peaklist_final))
          for(mi in 1:length(minidist)){
            minidist[mi] = SAX_mindist(Peaklist_final$Sequence[mi],Peaks$Sequence[t],Peaklist_final$Nr_of_Points[mi],SAX_reference_table)
          }
          Peaklist_final = Peaklist_final[match(min(minidist),minidist),]
          sample_peak = Peaklist_final
        } else {
          Peaklist_final = Peaklist_final[match(min(diff),diff),]
          sample_peak = Peaklist_final
        }
      }
    } else {
      shift_corrected_times = Peaks$RT[t]+Peaklist_final$intensity_shift+IS_dependent_shift_compound
      diff = abs(Peaklist_final$RT-shift_corrected_times)
      diff.grid = expand.grid(diff,diff)
      diff.grid$diff = abs(diff.grid$Var1-diff.grid$Var2)
      diff.grid = diff.grid[diff.grid$diff!=0,]
      differences = unique(diff.grid$diff)
      nearest_peaks = min(diff)==diff.grid$Var1[diff.grid$diff==min(diff.grid$diff)]
      if(min(differences)<0.05&T%in%nearest_peaks&use.MINDIST==T){
        minidist = rep(NA,nrow(Peaklist_final))
        for(mi in 1:length(minidist)){
          minidist[mi] = SAX_mindist(Peaklist_final$Sequence[mi],Peaks$Sequence[t],Peaklist_final$Nr_of_Points[mi],SAX_reference_table)
        }
        Peaklist_final = Peaklist_final[match(min(minidist),minidist),]
        sample_peak = Peaklist_final
      } else {
        Peaklist_final = Peaklist_final[match(min(diff),diff),]
        sample_peak = Peaklist_final
      }
    }
  } else {
    sample_peak = Peaklist_final
  }
  
  sample_peak$LOD = LOD
  sample_peak$LOQ = LOQ
  
  intensity_dependent_shift = 0
  if(Peaks$Compound[t]%in%Intensity_dependent_shift_names&nrow(sample_peak[!is.na(sample_peak$RT),])>0){
    frame_intensity_shift = Intensity_dependent_shift[[match(Peaks$Compound[t],Intensity_dependent_shift_names)]]
    frame_intensity_shift$diff = abs(sample_peak$Height-frame_intensity_shift$int)
    frame_intensity_shift$diff[frame_intensity_shift$diff==0] = 10^-100
    frame_intensity_shift$diff = 1/frame_intensity_shift$diff
    two_smallest = frame_intensity_shift[order(frame_intensity_shift$diff,decreasing = T),]
    total = sum(two_smallest$diff)
    intensity_dependent_shift = sum(two_smallest$shift*two_smallest$diff)/total
  }
  
  ID = substr(Peaks$Compound[t],1,4)
  if(T%in%grepl("IS",ID)){
    intensity_dependent_shift = 0
  }
  
  shift_corrected_time = quan_peak$RT+intensity_dependent_shift+IS_dependent_shift_compound
  shift_corrected_Start_level = quan_peak$Start_RT_level+intensity_dependent_shift+IS_dependent_shift_compound
  shift_corrected_End_level = quan_peak$End_RT_level+intensity_dependent_shift+IS_dependent_shift_compound
  
  if(bad_chroma==F&nrow(sample_peak[!is.na(sample_peak$Height),])>0){
    plot(EIC$int~EIC$rt,type="p",main=paste0(quan_peak$Compound))
    lines(EIC_smooth$rt,EIC_smooth$int,col="red",type = "p")
    for(qlk in 1:nrow(sample_peak)){
      abline(v=sample_peak$Start_RT[qlk],col="purple")
      abline(v=sample_peak$End_RT[qlk],col="purple")
      Peak = EIC[EIC$rt>=sample_peak$Start_RT[qlk]&EIC$rt<=sample_peak$End_RT[qlk],]
      lines(Peak$rt,Peak$int,type="l",col="purple",lwd=2)
    }
  }
  
  #Check if sequence is within acceptable MINDIST range
  mindist_too_large = F
  if(use.MINDIST==T){
    if(length(sample_peak$RT[!is.na(sample_peak$RT)])>0){
      mindist = plyr::round_any(SAX_mindist(sample_peak$Sequence,Peaks$Sequence[t],sample_peak$Nr_of_Points,SAX_reference_table),f=floor,accuracy = 0.1)
      maximum_mindist = plyr::round_any(Peaks$max_MINDIST[t]/sqrt(Peaks$Nr_of_Points[t])*sqrt(sample_peak$Nr_of_Points),f=ceiling,accuracy=0.1)
      
      if(mindist > maximum_mindist){
        mindist_too_large=T
        sample_peak$Comment = paste0(sample_peak$Comment," Peakshape differs too much_<a7>",sample_peak$Area,"<a7>;")
      }
    }
  } 
  Peak_sub = EIC[EIC$rt>=sample_peak$Start_RT&EIC$rt<=sample_peak$End_RT,]
  
  #Get criteria
  mindist_too_large = mindist_too_large
  bad_chroma = bad_chroma
  really_above_LOD = as.logical(((nrow(Peak_sub[Peak_sub$int>=LOD,])/nrow(Peak_sub))>minimum_qualitative_threshold)&((sample_peak$Height/quantile(EIC$int[EIC$rt<=sample_peak$Start_RT|EIC$rt>=sample_peak$End_RT],probs=c(0.95),na.rm=T))>minimum_background_ratio))
  above_LOQ = sample_peak$Height>LOQ
  Solvent_blank = grepl(Solvent_Blank_ID,files[j])
  if(allowed_shift>Peaks$RT_tol[t]){
    tolerance = allowed_shift
  } else {
    tolerance = Peaks$RT_tol[t]
  }
  #You have to allow all peaks to shift around the allowed shift, otherwise if the shift was calculated wrong (i.e. IS shifts, but analyte not) you would virtually add errors
  within_time = !(abs(sample_peak$RT-shift_corrected_time)>(allowed_shift+Peaks$RT_tol[t])|abs(sample_peak$Start_RT_level-shift_corrected_Start_level)>(Peaks$Start_var[t]+allowed_shift)|abs(sample_peak$End_RT_level-shift_corrected_End_level)>(Peaks$End_var[t]+allowed_shift))
  if(is.na(within_time)){
    within_time = F
  }
  within_CHECK_range = !(abs(sample_peak$RT-shift_corrected_time)>(allowed_shift+Peaks$RT_tol[t])|abs(sample_peak$Start_RT_level-shift_corrected_Start_level)>(Peaks$Start_var[t]+allowed_shift)|abs(sample_peak$End_RT_level-shift_corrected_End_level)>(Peaks$End_var[t]+allowed_shift))
  
  if(nrow(Peaks_list[[t]][Peaks_list[[t]]$type!="quan",])>0){
    fully_confirmed = sample_peak$confirmed==nrow(Peaks_list[[t]][Peaks_list[[t]]$type!="quan",])
    plausible_peak_wo_confirming = within_time&(above_LOQ|really_above_LOD)&sample_peak$confirmed!=nrow(Peaks_list[[t]][Peaks_list[[t]]$type!="quan",])&sample_peak$confirmed>=0
  } else {
    fully_confirmed = T
    plausible_peak_wo_confirming = F
  }
  
  #Check if a double peak emerged from the sample and the EIC needs to be wider
  if(within_time&!above_LOQ&sample_peak$Height/min(EIC$int[EIC$int!=0])>10&sample_search_window!=(width*width_factor_background)){
    left_baseline = median(EIC$int[EIC$rt<=quantile(EIC$rt,probs=c(0.1))],na.rm = T)<=Baseline*extended_baseline_factor
    right_baseline = median(EIC$int[EIC$rt>=quantile(EIC$rt,probs=c(0.9))],na.rm = T)<=Baseline*extended_baseline_factor
    if(!(left_baseline&right_baseline)&(left_baseline|right_baseline)){
      sample_analysis(sample_path,results_path,sample,files,Peaks,Peaks_list,Solvent_Blank_ID,IS_dependent_shift,Intensity_dependent_shift,Multiple,FullscanMS1,precursor_window,inner_intensity_ratio_multiple_peaks,Scan_filters,IS_Assignment,compounds_shift_corrected,RT_range,ppm_val,method_time,IS_deviation,gen.plots,use.MINDIST,SAX_reference_table,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,minimum_background_ratio,extended_baseline_factor,width_smoothing,width_factor_background=width_factor_background*2,sample_search_window,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,maximum_allowed_shift,maximum_allowed_shift_ratio,minimum_cutoff_intensity_factor,minimum_confirming_peak_height,max_MINDIST,minimum_datapoints_per_sample_peak,minimum_qualitative_threshold,t,j,rep=T)
    }
  }
  
  #Check if the shift of the Peak is too high compared to standards
  if(nrow(sample_peak[!is.na(sample_peak$RT),])>0){
    if((within_time&!mindist_too_large&!above_LOQ&really_above_LOD&(fully_confirmed|plausible_peak_wo_confirming))){
      sample_peak$Comment = paste0(sample_peak$Comment," only_above_LOD_?",sample_peak$Area,"?;")
      if(gen.plots==T){
        dir.create(paste0(results_path,"/plots"))
        if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
          Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
          PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        } else {
          PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        }
        pdf(file=PDF_path)
        p <- ggplot(EIC, aes(rt, int)) + 
          geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
        p = p + 
          labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
               subtitle = "ONLY_QUAL") +
          theme_classic() +
          theme(plot.title = element_text(lineheight = 1.1),
                plot.subtitle = element_text(lineheight = 1.0,colour = "lightgreen"),
                legend.position = "none")
        p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
        p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
        p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
        p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                           size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                         size=1.5)
        Peak = EIC[EIC$rt>=sample_peak$Start_RT&EIC$rt<=sample_peak$End_RT,]
        p=p + geom_line(data=Peak,aes(rt,int,col="purple"),size=2)
        p=p + geom_vline(xintercept=sample_peak$Start_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$End_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$RT,col="purple",
                         size=1.5)
        p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                                 labels = scales::number_format(accuracy = 0.01))
        print(p)
        dev.off()
      }
    } 
    else if((!within_time&within_CHECK_range&really_above_LOD&!Solvent_blank)|(within_time&sample_peak$Nr_of_Points<minimum_datapoints_per_sample_peak&above_LOQ&!Solvent_blank)){
      if((!within_time&within_CHECK_range&really_above_LOD&!Solvent_blank)){
        sample_peak$Comment = paste0(sample_peak$Comment," shift reference:sample was too high_<a7>",sample_peak$Area,"<a7>;")
        title = "CHECK - shift reference:sample was too high"
      } else if((within_time&sample_peak$Nr_of_Points<minimum_datapoints_per_sample_peak&above_LOQ&!Solvent_blank)){
        #few data points only if the EIC is rather clean
        plausibility_ratio = sample_peak$Height/median(EIC$int[EIC$rt<=sample_peak$Start_RT|EIC$rt>=sample_peak$End_RT],na.rm = T)
        
        if(plausibility_ratio>minimum_background_ratio){
          sample_peak$Comment = paste0(sample_peak$Comment," few data points_<a7>",sample_peak$Area,"<a7>;")
          title = "CHECK - only few data points"
        } else {
          title = "NOT FOUND"
        }
      } else if(plausible_peak_wo_confirming){
        sample_peak$Comment = paste0(sample_peak$Comment," confirming ion missing_<a7>",sample_peak$Area,"<a7>;")
        title = "CHECK - expected confirming ion(s) missing"
      } else {
        title = "CHECK"
      }
      if(gen.plots==T){
        dir.create(paste0(results_path,"/plots"))
        if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
          Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
          PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        } else {
          PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        }
        pdf(file=PDF_path)
        p <- ggplot(EIC, aes(rt, int)) + 
          geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
        p = p + 
          labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
               subtitle = title) +
          theme_classic() +
          theme(plot.title = element_text(lineheight = 1.1),
                plot.subtitle = element_text(lineheight = 1.0,colour = "orange"),
                legend.position = "none")
        p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
        p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
        p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
        p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                           size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                         size=1.5)
        Peak = EIC[EIC$rt>=sample_peak$Start_RT&EIC$rt<=sample_peak$End_RT,]
        p=p + geom_line(data=Peak,aes(rt,int,col="purple"),size=2)
        p=p + geom_vline(xintercept=sample_peak$Start_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$End_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$RT,col="purple",
                         size=1.5)
        p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                                 labels = scales::number_format(accuracy = 0.01))
        print(p)
        dev.off()
      }
      sample_peak[1,4:(length(sample_peak)-2)] = NA
    } 
    else if((within_time&!mindist_too_large&above_LOQ)|(plausible_peak_wo_confirming)){
      sample_peak = sample_peak
      if(gen.plots==T){
        dir.create(paste0(results_path,"/plots"))
        if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
          Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
          PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        } else {
          PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        }
        if((plausible_peak_wo_confirming)){
          title = "FOUND_BUT_NOT_CONFIRMED"
        } else {
          title = "FOUND_AND_CONFIRMED"
        }
        pdf(file=PDF_path)
        p <- ggplot(EIC, aes(rt, int)) + 
          geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
        p = p + 
          labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
               subtitle = title) +
          theme_classic() +
          theme(plot.title = element_text(lineheight = 1.1),
                plot.subtitle = element_text(lineheight = 1.0,colour = "darkgreen"),
                legend.position = "none")
        p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
        p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
        p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
        p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                           size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                         size=1.5)
        Peak = EIC[EIC$rt>=sample_peak$Start_RT&EIC$rt<=sample_peak$End_RT,]
        p=p + geom_line(data=Peak,aes(rt,int,col="purple"),size=2)
        p=p + geom_vline(xintercept=sample_peak$Start_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$End_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$RT,col="purple",
                         size=1.5)
        p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                                 labels = scales::number_format(accuracy = 0.01))
        print(p)
        dev.off()
      }
    } 
    else {
      if(grepl("confirming",sample_peak$Comment)){
        if(gen.plots==T){
          dir.create(paste0(results_path,"/plots"))
          if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
            Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
            PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
          } else {
            PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
          }
          pdf(file=PDF_path)
          p <- ggplot(EIC, aes(rt, int)) + 
            geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
          p = p + 
            labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
                 subtitle = "CONFIRMING_MISSING") +
            theme_classic() +
            theme(plot.title = element_text(lineheight = 1.1),
                  plot.subtitle = element_text(lineheight = 1.0,colour = "brown"),
                  legend.position = "none")
          p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
          p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
          p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
          p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                             size=1.5)
          p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                           size=1.5)
          p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                           size=1.5)
          p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                           size=1.5)
          p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                                   labels = scales::number_format(accuracy = 0.01))
          print(p)
          dev.off()
        }
      }
      if(gen.plots==T){
        dir.create(paste0(results_path,"/plots"))
        if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
          Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
          PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        } else {
          PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
        }
        pdf(file=PDF_path)
        p <- ggplot(EIC, aes(rt, int)) + 
          geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
        p = p + 
          labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
               subtitle = "NOT_FOUND") +
          theme_classic() +
          theme(plot.title = element_text(lineheight = 1.1),
                plot.subtitle = element_text(lineheight = 1.0,colour = "darkgrey"),
                legend.position = "none")
        p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
        p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
        p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
        p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                           size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                         size=1.5)
        Peak = EIC[EIC$rt>=sample_peak$Start_RT&EIC$rt<=sample_peak$End_RT,]
        p=p + geom_line(data=Peak,aes(rt,int,col="purple"),size=2)
        p=p + geom_vline(xintercept=sample_peak$Start_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$End_RT_level,col=adjustcolor(col="purple",alpha.f=0.5),
                         size=1.5)
        p=p + geom_vline(xintercept=sample_peak$RT,col="purple",
                         size=1.5)
        p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                                 labels = scales::number_format(accuracy = 0.01))
        print(p)
        dev.off()
      }
      sample_peak[1,4:(length(sample_peak)-2)] = NA
    }
  } else if(bad_chroma){
    
    if (T%in%grepl(".mzML", files)) {
      if(MS_filter==FullscanMS1){
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "EIC",
                               mz = quan_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
      } else {
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "MS2",
                               mz = quan_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        precursor = quan_peak$mz
        EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
        EIC = EIC[(EIC$mz>=RaMS::pmppm(quan_peak$mz,ppm_val)[1]&EIC$mz<=RaMS::pmppm(quan_peak$mz,ppm_val)[2]),]
        EIC = EIC[,c(1:2)]
      }
    } else {
      
      Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                       mass = quan_peak$mz,
                                       tol = ppm_val*2,
                                       filter = MS_filter,
                                       type = "xic")
      
      EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
    }
    
    EIC_background = EIC
    EIC = EIC[EIC$rt>=(shift_corrected_Start_level-5*RT_range)&EIC$rt<=(shift_corrected_End_level+5*RT_range),]
    EIC = EIC[!is.na(EIC$rt),]
    if(grepl(".mzML",files)){
      if(nrow(EIC)<3){
        EIC = data.frame("rt"=c(0,method_time/3,method_time*2/3,method_time),"int"=c(0,0,0,0))
      }
    } else {
      if(nrow(EIC)<3){
        EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
      }
    }
    timesplit = ceiling(length(EIC$rt)/3)
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(i in 1:3){
      EICs[[i]] = EIC[(Start:(Start+timesplit)),]
      EICs[[i]] = EICs[[i]][!is.na(EICs[[i]]$rt)&!is.na(EICs[[i]]$int),]
      timerange = max(EICs[[i]]$rt)-min(EICs[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_smooth[[i]] = as.data.frame(ksmooth(EICs[[i]]$rt,EICs[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[i]] = EICs_smooth[[i]][!is.na(EICs_smooth[[i]]$x)&!is.na(EICs_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    
    sample_peak$Comment = paste0(sample_peak$Comment," not_found_&",sample_peak$Area,"&;")
    if(gen.plots==T){
      dir.create(paste0(results_path,"/plots"))
      if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
        Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
        PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      } else {
        PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      }
      pdf(file=PDF_path)
      p <- ggplot(EIC, aes(rt, int)) + 
        geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
      p = p + 
        labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
             subtitle = "BAD_CHROMATOGRAM_NOT_FOUND") +
        theme_classic() +
        theme(plot.title = element_text(lineheight = 1.1),
              plot.subtitle = element_text(lineheight = 1.0,colour = "darkgrey"),
              legend.position = "none")
      p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                       size=1.5)
      p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                       size=1.5)
      p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                       size=1.5)
      p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                               labels = scales::number_format(accuracy = 0.01))
      print(p)
      dev.off()
    }
  } else {
    
    if (T%in%grepl(".mzML", files)) {
      if(MS_filter==FullscanMS1){
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "EIC",
                               mz = quan_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        EIC = data.frame("rt"=EIC$EIC$rt,"int"=EIC$EIC$int)
      } else {
        EIC = RaMS::grabMSdata(files = paste0(sample_path,"/",files[j]),
                               grab_what = "MS2",
                               mz = quan_peak$mz,
                               ppm = ppm_val,
                               rtrange = NULL,
                               prefilter = -1)
        precursor = quan_peak$mz
        EIC = EIC[((EIC$prec>=(precursor-precursor_window))&(EIC$prec<=(precursor+precursor_window))),]
        EIC = EIC[(EIC$mz>=RaMS::pmppm(quan_peak$mz,ppm_val)[1]&EIC$mz<=RaMS::pmppm(quan_peak$mz,ppm_val)[2]),]
        EIC = EIC[,c(1:2)]
      }
    } else {
      
      Chroma = rawrr::readChromatogram(rawfile = paste0(sample_path,"/",files[j]),
                                       mass = quan_peak$mz,
                                       tol = ppm_val*2,
                                       filter = MS_filter,
                                       type = "xic")
      
      EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
    }
    EIC_background = EIC
    EIC = EIC[EIC$rt>=(shift_corrected_Start_level-1)&EIC$rt<=(shift_corrected_End_level+1),]
    EIC = EIC[!is.na(EIC$rt),]
    if(grepl(".mzML",files)){
      if(nrow(EIC)<3){
        EIC = data.frame("rt"=c(0,method_time/3,method_time*2/3,method_time),"int"=c(0,0,0,0))
      }
    } else {
      if(nrow(EIC)<3){
        EIC = data.frame("rt"=Chroma[[1]]$times,"int"=Chroma[[1]]$intensities)
      }
    }
    
    timesplit = ceiling(length(EIC$rt)/3)
    EICs = list()
    EICs_smooth = list()
    Start = 1
    for(i in 1:3){
      EICs[[i]] = EIC[(Start:(Start+timesplit)),]
      EICs[[i]] = EICs[[i]][!is.na(EICs[[i]]$rt)&!is.na(EICs[[i]]$int),]
      timerange = max(EICs[[i]]$rt)-min(EICs[[i]]$rt)
      bandwidth = 0.1*timerange
      EICs_smooth[[i]] = as.data.frame(ksmooth(EICs[[i]]$rt,EICs[[i]]$int,kernel="normal",bandwidth = bandwidth))
      EICs_smooth[[i]] = EICs_smooth[[i]][!is.na(EICs_smooth[[i]]$x)&!is.na(EICs_smooth[[i]]$y),]
      Start=Start+timesplit
      if(Start>length(EIC$rt)){
        break
      }
    }
    EIC_smooth = dplyr::bind_rows(EICs_smooth)
    colnames(EIC_smooth) = c("rt","int")
    
    sample_peak$Comment = paste0(sample_peak$Comment," not_found_&",sample_peak$Area,"&;")
    if(gen.plots==T){
      dir.create(paste0(results_path,"/plots"))
      if(nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))>260){
        Diff_PDF = 250 - nchar(paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf"))
        PDF_path = paste0(results_path,"/plots/",stringr::str_sub(stringr::str_remove(files[j],".raw"),1,nchar(stringr::str_remove(files[j],".raw")+Diff_PDF)),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      } else {
        PDF_path = paste0(results_path,"/plots/",stringr::str_remove(files[j],".raw"),"___",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[t])),".pdf")
      }
      pdf(file=PDF_path)
      p <- ggplot(EIC, aes(rt, int)) + 
        geom_point(size = 2) + xlab("Retention Time [Min]") + ylab("Intensity")
      p = p + 
        labs(title = stringr::str_wrap(paste0(files[j]," - ",quan_peak$Compound),width=70), 
             subtitle = "NOT_FOUND") +
        theme_classic() +
        theme(plot.title = element_text(lineheight = 1.1),
              plot.subtitle = element_text(lineheight = 1.0,colour = "darkgrey"),
              legend.position = "none")
      p=p + geom_point(data = EIC_smooth,aes(rt,int,col="red"),size=2)
      p = p + geom_hline(yintercept=LOD,col="lightgreen",size=1.5)
      p = p + geom_hline(yintercept=LOQ,col="darkgreen",size=1.5)
      p = p + geom_hline(yintercept=Baseline,col=adjustcolor(col="darkgrey",alpha.f = 0.8),
                         size=1.5)
      p=p + geom_vline(xintercept=shift_corrected_Start_level,col=adjustcolor(col="black",alpha.f = 0.5),
                       size=1.5)
      p=p + geom_vline(xintercept=shift_corrected_End_level,col=adjustcolor(col="black",alpha.f = 0.5),
                       size=1.5)
      p=p + geom_vline(xintercept=shift_corrected_time,col="black",
                       size=1.5)
      p = p+scale_x_continuous(breaks=seq(min(EIC$rt),max(EIC$rt),by=(max(EIC$rt)-min(EIC$rt))/10),
                               labels = scales::number_format(accuracy = 0.01))
      print(p)
      dev.off()
    }
  }
  
  return(sample_peak)
}

#Compiled functions
################################################################################
#'setup_environment
#'
#'sets up the working environment within the R session
#'
#' @param sample the name of the sample / file batch
#' @export
setup_environment = function(sample){
  currentFunction = "setup of environment"
  sample_paths = normalizePath(gtools::mixedsort(list.files("./2_Samples",full.names = T)[]))
  s = match(sample,gtools::mixedsort(list.files("./2_Samples")))
  sample_path = sample_paths[s]
  results_path = paste0("./3_Results/",sample)
  if(file.exists(stringr::str_replace(sample_path,"2_Samples","3_Results"))){
    if(length(list.files(stringr::str_replace(sample_path,"2_Samples","3_Results"),pattern = "Results_"))==0){
      print(paste0(sample," ready for analysis"))
      analyze = T
    } else {
      print(paste0(sample," skipped"))
      analyze = F
    }
  } else {
    print(paste0(sample," ready for analysis"))
    analyze = T
  }
  
  #Settings from YAML file
  settings = yaml::yaml.load_file(paste0(list.files(sample_path,pattern="yaml",full.names = T)))
  list2env(settings,envir = .GlobalEnv)
  
  #Load target files
  compounds = data.table::fread(paste0("./1_Targets/",Targets_file))
  if(!is.null(Qual_ID)){
    compounds_qual = data.table::fread(paste0("./1_Targets/",Qual_file))
  } else {
    compounds_qual = NA
  }
  IS = compounds[grepl("IS",compounds$ID),] #used for IS-dependent analysis
  IS = IS[!grepl(paste0(letters,collapse = "|"),IS$ID),]
  IS_Assignment = data.table::fread(paste0("./1_Targets/",IS_Assignment_file))
  Multiple = data.table::fread(paste0("./1_Targets/",Multiple_Peaks_file))
  
  #.raw or .mzml files to be analyzed
  files = list.files(sample_path)
  files = files[grepl(".raw|.mzML",files)]
  files = mixedsort(files)
  
  #Define References
  References = files[grepl(Reference_ID,files)]
  
  #Create SAX reference table
  SAX_reference_table = jmotif::sax_distance_matrix(alphabet_size+1)
  colnames(SAX_reference_table) = letters[1:(alphabet_size+1)]
  rownames(SAX_reference_table) = colnames(SAX_reference_table)
  
  settings$IS = IS
  settings$files = files
  settings$References = References
  settings$sample_path = sample_path
  settings$results_path = results_path
  settings$compounds = compounds
  settings$compounds_qual = compounds_qual
  settings$IS_Assignment = IS_Assignment
  settings$Multiple = Multiple
  settings$analyze = analyze
  settings$SAX_reference_table = SAX_reference_table
  return(settings)
}

#'calculate_system_dependent_shift_and_RT_range
#'
#'Calculates the system-dependent shift and defined the retention time range
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @param minimum_search_window the smallest size of the search window in minutes
#' @export
calculate_system_dependent_shift_and_RT_range = function(sample,set.env=F,minimum_search_window){
  if(set.env==T){
    environment(setup_environment)= .GlobalEnv
    list2env(setup_environment(sample),envir=.GlobalEnv)
  }
  
  pb <- txtProgressBar(max=length(References), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  cl <- makeCluster(length(References))
  doSNOW::registerDoSNOW(cl)
  currentFunction = "calculating system-dependent shift"
  environment(errors) = environment()
  tryCatch({
    print("Calculating system-dependent shift")
    shift_all = foreach(q = 1:length(References),
                        .packages = c("seewave","dplyr","rawrr","foreach","doSNOW","zoo","parallel","RaMS"),
                        .options.snow=opts,
                        .combine = "c",
                        .export=c("References","sample_path","IS","FullscanMS1","ppm_val","cores","functions",names(functions),"method_time","alphabet_size","zigzag_trigger_threshold","extended_baseline_factor","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","minimum_peak_area","maximum_peak_width")) %dopar% {
                          
                          assign(paste0("cl",q,collapse = ""),makeCluster((cores-length(References))/length(References)))
                          
                          doSNOW::registerDoSNOW(get(paste0("cl",q,collapse = "")))
                          
                          shift_all = foreach(i = 1:length(IS$ID),
                                              .packages = c("seewave","dplyr","zoo","rawrr","RaMS"),
                                              .options.snow=opts,
                                              .combine = "c",
                                              .export=c("References","sample_path","IS","FullscanMS1","ppm_val","functions",names(functions),"q","method_time","alphabet_size","zigzag_trigger_threshold","extended_baseline_factor","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","minimum_peak_area","maximum_peak_width")) %dopar% {
                                                calc_shift(sample_path,References,IS,ppm_val,method_time,FullscanMS1,alphabet_size,zigzag_trigger_threshold,extended_baseline_factor,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,minimum_peak_area,maximum_peak_width,i,q)
                                              }
                          
                          stopCluster(get(paste0("cl",q,collapse = "")))
                          
                          shift_all
                        }
    close(pb)
    stopCluster(cl)
    
    shift_all = shift_all[!is.na(shift_all)]
    shift_all = shift_all[!abs(shift_all)>2]
    
    Shift = median(shift_all)
    
    RT_range = abs(max(shift_all)-min(shift_all))
    RT_range_rounded = round(RT_range,digits=1)
    if(RT_range_rounded<RT_range){
      RT_range = RT_range_rounded+0.1
    } else {
      RT_range = RT_range_rounded
    }
    
    if(T%in%grepl("pos|neg",sample)){
      if(RT_range<0.5){
        RT_range = 0.5
      }
    } else {
      if(RT_range<0.1){
        RT_range = 0.1
      }
    }
    
    RT_range = as.numeric(RT_range)
    
    print(paste0("system-dependent shift is ",round(Shift,digits = 2)," min and search window is ",round(RT_range,digits = 2)," min"))
    
    compounds_shift_corrected = compounds
    compounds_shift_corrected$`retention time` = compounds$`retention time`+Shift
    compounds_without_qual = compounds_shift_corrected[!grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID),]
    
    if(!is.null(Qual_ID)){
      compounds_qual_shift_corrected = compounds_qual
      compounds_qual_shift_corrected$`retention time` = compounds_qual$`retention time`+Shift
      compounds_qual_without_qual = compounds_qual_shift_corrected[!grepl(paste0(letters,collapse = "|"),compounds_qual_shift_corrected$ID),]
    }
    
    Output = list()
    Output$Shift = Shift
    Output$RT_range = RT_range
    Output$compounds_shift_corrected = compounds_shift_corrected
    Output$compounds_without_qual = compounds_without_qual
    if(!is.null(Qual_ID)){
      Output$compounds_qual_shift_corrected = compounds_qual_shift_corrected
      Output$compounds_qual_without_qual = compounds_qual_without_qual
    }
    return(Output)
  },error = errors)
  
}

#'calculate_system_dependent_shift_and_RT_range_woParallel
#'
#'Calculates the system-dependent shift and defines the retention time range but without parallelization
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @param minimum_search_window the smallest size of the search window in minutes
#' @export
calculate_system_dependent_shift_and_RT_range_woParallel = function(sample,set.env=F,minimum_search_window){
  if(set.env==T){
    environment(setup_environment) = .GlobalEnv
    list2env(setup_environment(sample),envir=.GlobalEnv)
  }
  
  currentFunction = "calculating system-dependent shift"
  pb <- txtProgressBar(max=nrow(IS), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  q = 1
  
  list2env(functions)
  
  cl <- makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  print("Calculating system-dependent shift")
  environment(errors) = environment()
  
  tryCatch({
    print(paste0("Checking: ",References[q]))
    shift_all = foreach::foreach(i = 1:length(IS$ID),
                                 .packages = c("seewave","dplyr","doSNOW","parallel","zoo","rawrr","stats","RaMS"),
                                 .options.snow=opts,
                                 .combine = "c",
                                 .export=c("References","sample_path","IS","FullscanMS1","ppm_val","functions",names(functions),"q","method_time","alphabet_size","zigzag_trigger_threshold","extended_baseline_factor","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","minimum_peak_area","maximum_peak_width")) %dopar% {
                                   calc_shift(sample_path,References,IS,ppm_val,method_time,FullscanMS1,alphabet_size,zigzag_trigger_threshold,extended_baseline_factor,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,minimum_peak_area,maximum_peak_width,i,q)
                                 }
  },error = errors)
  
  close(pb)
  stopCluster(cl)
  
  if(length(References)>1){
    for(q in 2:length(References)){
      print(paste0("Checking: ",References[q]))
      pb <- txtProgressBar(max=nrow(IS), style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      
      cl <- makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
      environment(errors) = environment()
      tryCatch({
        shift_all_q = foreach::foreach(i = 1:length(IS$ID),
                                       .packages = c("seewave","dplyr","doSNOW","parallel","zoo","rawrr","RaMS"),
                                       .options.snow=opts,
                                       .combine = "c",
                                       .export=c("References","sample_path","IS","FullscanMS1","ppm_val","functions",names(functions),"q","method_time","alphabet_size","zigzag_trigger_threshold","extended_baseline_factor","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","minimum_peak_area","maximum_peak_width")) %dopar% {
                                         calc_shift(sample_path,References,IS,ppm_val,method_time,FullscanMS1,alphabet_size,zigzag_trigger_threshold,extended_baseline_factor,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,minimum_peak_area,maximum_peak_width,i,q)
                                       }
      },error=errors)
      
      close(pb)
      stopCluster(cl)
      
      shift_all = c(shift_all,shift_all_q)
    }
  }
  
  shift_all = shift_all[!is.na(shift_all)]
  shift_all = shift_all[!abs(shift_all)>2]
  
  Shift = median(shift_all)
  
  RT_range = abs(max(shift_all)-min(shift_all))
  RT_range_rounded = round(RT_range,digits=1)
  if(RT_range_rounded<RT_range){
    RT_range = RT_range_rounded+0.1
  } else {
    RT_range = RT_range_rounded
  }
  
  if(T%in%grepl("pos|neg",sample)){
    if(RT_range<0.5){
      RT_range = 0.5
    }
  } else {
    if(RT_range<0.1){
      RT_range = 0.1
    }
  }
  
  RT_range = as.numeric(RT_range)
  
  print(paste0("system-dependent shift is ",round(Shift,digits = 2)," min and search window is ",round(RT_range,digits = 2)," min"))
  
  compounds_shift_corrected = compounds
  compounds_shift_corrected$`retention time` = compounds$`retention time`+Shift
  compounds_without_qual = compounds_shift_corrected[!grepl(paste0(letters,collapse = "|"),compounds_shift_corrected$ID),]
  
  if(!is.null(Qual_ID)){
    compounds_qual_shift_corrected = compounds_qual
    compounds_qual_shift_corrected$`retention time` = compounds_qual$`retention time`+Shift
    compounds_qual_without_qual = compounds_qual_shift_corrected[!grepl(paste0(letters,collapse = "|"),compounds_qual_shift_corrected$ID),]
  }
  
  Output = list()
  Output$Shift = Shift
  Output$RT_range = RT_range
  Output$compounds_shift_corrected = compounds_shift_corrected
  Output$compounds_without_qual = compounds_without_qual
  if(!is.null(Qual_ID)){
    Output$compounds_qual_shift_corrected = compounds_qual_shift_corrected
    Output$compounds_qual_without_qual = compounds_qual_without_qual
  }
  return(Output)
}

#'generate_peaklist_and_calculate_intensity_dependent_shift
#'
#'Generates the peaklist derived from references and checks for intensity-dependent shifts in calibrations
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @param use_shifting_search_window boolean, should the search window extend if intensities are detected at the edges of the EIC?
#' @param minimum_RT_tolerance the minimum tolerance in minutes of apex retention times of a sample peak compared to expectation
#' @param minimum_RT_tolerance_start the minimum tolerance in minutes of starting retention times of a sample peak compared to expectation
#' @param minimum_RT_tolerance_end the minimum tolerance in minutes of ending retention times of a sample peak compared to expectation
#' @export
generate_peaklist_and_calculate_intensity_dependent_shift = function(sample,set.env=F,use_shifting_search_window,maximum_RT_tolerance,maximum_RT_tolerance_start,maximum_RT_tolerance_end){
  
  if(set.env==T){
    environment(calculate_system_dependent_shift_and_RT_range) = .GlobalEnv
    environment(calculate_system_dependent_shift_and_RT_range_woParallel) = .GlobalEnv
    if(cores<2*length(files[grepl(Reference_ID,files)])){
      list2env(calculate_system_dependent_shift_and_RT_range_woParallel(sample,set.env=F),envir=.GlobalEnv)
    } else {
      list2env(calculate_system_dependent_shift_and_RT_range(sample,set.env=F),envir=.GlobalEnv)
    }
  }
  
  Peaklists = list()
  
  References_refined = gsub(paste0(".*",Reference_ID),Reference_ID,References)
  References_refined = substr(References_refined,(nchar(Reference_ID)+2),(nchar(Reference_ID)+5))
  References_refined = gsub("p",".",References_refined)
  References_sorted = References_refined[order(References_refined,decreasing = T)]
  References_refined = References[(grepl(paste0(References_sorted,collapse = "|"),References))]
  
  Peaks_all_list = list()
  
  cl <- makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  
  currentFunction = "generate peaklist"
  
  environment(errors) = environment()
  tryCatch({
    for(Reference in 1:length(References_refined)){
      
      print(paste0("Generate peaklist of reference: ",References_refined[Reference]))
      
      pb <- txtProgressBar(max=length(compounds_shift_corrected$ID), style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      
      Peaks_all = foreach(z=1:length(compounds_shift_corrected$ID),
                          .packages = c("seewave","dplyr","rawrr","zoo","stats","MESS","tcpl","RaMS"),
                          .combine = "rbind",
                          .multicombine = T,
                          .options.snow=opts,
                          .export = c("sample_path",
                                      "References_refined",
                                      "compounds_shift_corrected",
                                      "ppm_val",
                                      "Reference",
                                      "RT_range",
                                      names(functions),
                                      "FullscanMS1",
                                      "precursor_window",
                                      "outer_search_window_multiple_peaks",
                                      "inner_search_window_multiple_peaks",
                                      "outer_intensity_ratio_multiple_peaks",
                                      "inner_intensity_ratio_multiple_peaks",
                                      "Scan_filters",
                                      "Multiple",
                                      "sample",
                                      "method_time",
                                      "alphabet_size",
                                      "use_shifting_search_window",
                                      "zigzag_trigger_threshold",
                                      "extended_baseline_factor",
                                      "normal_background_quantile",
                                      "higher_background_quantile",
                                      "minimum_background_ratio",
                                      "width_smoothing",
                                      "maximum_nr_of_a",
                                      "minimum_nr_of_high_intensity_letters",
                                      "maximum_peak_width",
                                      "minimum_peak_area",
                                      "minimum_nr_of_datapoints_per_peak")) %dopar% {
                                        get_peaks_all(sample_path,sample,RT_range,Multiple,method_time,FullscanMS1,Scan_filters,precursor_window,outer_search_window_multiple_peaks,inner_search_window_multiple_peaks,outer_intensity_ratio_multiple_peaks,inner_intensity_ratio_multiple_peaks,References_refined,compounds_shift_corrected,ppm_val,Times_shift_r=0,Times_shift_l=0,use_shifting_search_window=F,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,extended_baseline_factor,minimum_background_ratio,width_smoothing,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,z,Reference)
                                      }
      
      Peaks_all_list[[Reference]] = Peaks_all
    }
  },error = errors)
  
  Peaks_all = bind_rows(Peaks_all_list)
  
  #If qualitative analytes are included -> Analyze qualitative_list and add to Peaks_all
  if(!is.null(Qual_ID)){
    Peaks_all_qual_list = list()
    
    Qual_files = files[grepl(Qual_ID,files)]
    
    currentFunction = "generate qualitative peaklist"
    
    environment(errors) = environment()
    tryCatch({
      for(Reference in 1:length(Qual_files)){
        print(paste0("Qualitative References analyzed: ",Qual_files[Reference]))
        pb <- txtProgressBar(max=length(compounds_qual_shift_corrected$ID), style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
        
        Peaks_all_qual = foreach(z=1:length(compounds_qual_shift_corrected$ID),
                                 .packages = c("seewave","dplyr","rawrr","zoo","stats","MESS","tcpl","RaMS"),
                                 .combine = "rbind",
                                 .multicombine = T,
                                 .options.snow=opts,
                                 .export = c("sample_path",
                                             "Qual_files",
                                             "compounds_qual_shift_corrected",
                                             "ppm_val",
                                             "Reference",
                                             "RT_range",
                                             names(functions),
                                             "FullscanMS1",
                                             "precursor_window",
                                             "outer_search_window_multiple_peaks",
                                             "inner_search_window_multiple_peaks",
                                             "outer_intensity_ratio_multiple_peaks",
                                             "inner_intensity_ratio_multiple_peaks",
                                             "Scan_filters",
                                             "Scan_filters",
                                             "Multiple",
                                             "sample",
                                             "method_time",
                                             "alphabet_size",
                                             "use_shifting_search_window",
                                             "zigzag_trigger_threshold",
                                             "extended_baseline_factor",
                                             "normal_background_quantile",
                                             "higher_background_quantile",
                                             "minimum_background_ratio",
                                             "width_smoothing",
                                             "maximum_nr_of_a",
                                             "minimum_nr_of_high_intensity_letters",
                                             "maximum_peak_width",
                                             "minimum_peak_area",
                                             "minimum_nr_of_datapoints_per_peak")) %dopar% {
                                               get_peaks_all(sample_path,sample,RT_range,Multiple,method_time,FullscanMS1,Scan_filters,precursor_window,outer_search_window_multiple_peaks,inner_search_window_multiple_peaks,outer_intensity_ratio_multiple_peaks,inner_intensity_ratio_multiple_peaks,References_refined=Qual_files,compounds_shift_corrected=compounds_qual_shift_corrected,ppm_val,Times_shift_r=0,Times_shift_l=0,use_shifting_search_window=use_shifting_search_window,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,extended_baseline_factor,minimum_background_ratio,width_smoothing,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,z,Reference)
                                             }
      }
    },error=errors)
    
    Peaks_all_qual_list[[Reference]] = Peaks_all_qual
    
    Peaks_all = rbind(Peaks_all,Peaks_all_qual)
    compounds_without_qual = rbind(compounds_without_qual,compounds_qual_without_qual)
    compounds_shift_corrected = rbind(compounds_shift_corrected,compounds_qual_shift_corrected)
  }
  
  close(pb)
  stopCluster(cl)
  
  Peaks = Peaks_all[1:length(compounds_without_qual$ID),]
  Peaks$Start_var = NA
  Peaks$End_var = NA
  Peaks$max_MINDIST = NA
  Peaks$Start_RT_level = NA
  Peaks$End_RT_level = NA
  
  Peaks$Compound = compounds_without_qual$identity
  Peaks$mz = compounds_without_qual$`m/z`
  Peaks[,3:length(Peaks)] = NA
  
  Peaks_list = list()
  
  currentFunction = "summarize peaklist"
  
  for(p in 1:length(Peaks$Compound)){
    zu = NA
    ID = substr(Peaks$Compound[p],1,4)
    selected = Peaks_all[grepl(ID,Peaks_all$Compound),]
    ions = as.character(as.data.frame(table(selected$Compound))$Var1)
    QuanIon = Peaks$Compound[p]
    QuanIon_peaks = Peaks_all[Peaks_all$Compound%in%QuanIon,]
    QuanIon_peaks = QuanIon_peaks[!is.na(QuanIon_peaks$RT),]
    QualIons = ions[ions!=QuanIon]
    QualIon_peaks = Peaks_all[Peaks_all$Compound%in%QualIons,]
    QualIon_peaks = QualIon_peaks[!is.na(QualIon_peaks$RT),]
    if(length(QualIons)>=2){
      QualIon_peaks = Peaks_all[Peaks_all$Compound==QualIons[1],]
      for(v in 2:length(QualIons)){
        QualIon_peaks = rbind(QualIon_peaks,Peaks_all[Peaks_all$Compound==QualIons[v],])
      }
      QualIon_peaks = QualIon_peaks[!is.na(QualIon_peaks$RT),]
    }
    if(length(QuanIon_peaks$Compound)==0){
      if(length(QualIon_peaks$Compound)!=0){
        if(length(QualIons)==1){
          Peaks[p,]$Comment = "quan not found, but qual"
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        } else {
          ionas = as.character(as.data.frame(table(QualIon_peaks$Compound))$Var1)
          if(length(ionas)!=length(QualIons)){
            Peaks[p,]$Comment = "neither quan nor all qual identified"
            final_ions = Peaks[p,]
            final_ions$ratio = 1
            final_ions$type = "quan"
            Peaks_list[[p]] = final_ions
            next
          } else {
            peaks_a = QualIon_peaks[QualIon_peaks$Compound==QualIons[1],]
            peaks_b = QualIon_peaks[QualIon_peaks$Compound==QualIons[2],]
            diff_a = 1:length(peaks_a$Compound)
            for(a in 1:length(peaks_a$Compound)){
              diff_a[a] = abs(peaks_a$RT[a]-peaks_b$RT)
            }
            diff_b = 1:length(peaks_b$Compound)
            for(b in 1:length(peaks_b$Compound)){
              diff_b[b] = abs(peaks_b$RT[b]-peaks_a$RT)
            }
            diff_ab = min(diff_a)
            if(diff_ab<=0.1){
              Int_a = peaks_a[match(max(peaks_a$Height),peaks_a$Height),]$Height
              Int_b = peaks_b[match(max(peaks_b$Height),peaks_b$Height),]$Height
              Peaks$Comment[p] = "quan not found, first two qual found, highest used"
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            } else {
              Peaks$Comment[p] = "first two qual found, but share no RT"
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
          }
        }
      } 
    }
    if(length(QuanIon_peaks$Compound)!=0){
      if(length(QualIons)==0){
        aimed_time = compounds_without_qual$`retention time`[p]
        diffi = abs(QuanIon_peaks$RT-aimed_time)
        Remaining = QuanIon_peaks[diffi<=RT_range,]
        Remaining_RT_1 = plyr::round_any(Remaining$RT,0.05,f=floor)
        if(grepl("GC",sample)){
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 2)))
        } else {
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 1)))
        }
        
        if(length(Remaining_RT$Var1)>1){
          Remaining_RT$intensities = 1:length(Remaining_RT$Var1)
          for(v in 1:length(Remaining_RT$Var1)){
            if(grepl("GC",sample)){
              Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[v]],na.rm = T)
            } else {
              Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[v]],na.rm = T)
            }
          }
          Remaining_RT = Remaining_RT[!is.nan(Remaining_RT$intensities),]
          if(grepl("GC",sample)){
            Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)-as.numeric(paste(Remaining_RT$Var1))),digits = 2)
          } else {
            Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)-as.numeric(paste(Remaining_RT$Var1))),digits = 1)
          }
          
          Remaining_RT$ratio_multiple_threshold = NA
          Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff<=inner_search_window_multiple_peaks] = 1/inner_intensity_ratio_multiple_peaks
          Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff>inner_search_window_multiple_peaks] = 1/outer_intensity_ratio_multiple_peaks
          
          Ratios = max(Remaining_RT$intensities)/Remaining_RT$intensities
          Rem_Ratios = Ratios[Ratios!=1]
          
          if(as.numeric(ID)%in%Multiple$ID){
            zu = match(as.numeric(ID),Multiple$ID)
            if(T%in%(Ratios>Remaining_RT$ratio_multiple_threshold)){
              Remaining_RT = Remaining_RT[Ratios<=Remaining_RT$ratio_multiple_threshold,]
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
              }
            }
            if(Multiple$multiple[zu]!=0){
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[Multiple$multiple[zu]],]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[Multiple$multiple[zu]],]
              }
            }
          } else if(mean(Rem_Ratios>=(1/outer_intensity_ratio_multiple_peaks),na.rm = T)==1){
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            }
          } else if((nrow(Remaining_RT[Remaining_RT$time_diff>=min(Remaining_RT$time_diff),])>1)){
            Remaining_RT = Remaining_RT[Remaining_RT$time_diff<=min(Remaining_RT$time_diff,na.rm=T),]
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            }
          } else {
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
            }
          }
        } 
        
        if(length(Remaining_RT$Var1)==0){
          Peaks$Comment[p] = paste0(Peaks$Comment[p]," RTs do not match;")
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        }
        
        if(length(Remaining$Compound)==0){
          Peaks$Comment[p] = paste0(Peaks$Comment[p]," quan ion not found;")
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        }
        if(!is.na(zu)){
          if(Multiple$multiple[zu]==0){
            Peaks[p,] = Remaining[1,]
            Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
            Peaks$Start_RT[p] = min(Remaining$Start_RT,na.rm = T)
            Peaks$Start_RT_level[p] = min(Remaining$Start_RT_level,na.rm = T)
            Peaks$End_RT_level[p] = max(Remaining$End_RT_level,na.rm=T)
            Peaks$RT[p] = median(Remaining$RT,na.rm = T)
            Peaks$End_RT[p] = max(Remaining$End_RT,na.rm = T)
            Peaks$Height[p] = max(Remaining$Height,na.rm = T)
            Peaks$Width[p] = Peaks$End_RT[p]-Peaks$Start_RT[p]
            Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
            Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
            Peaks$density[p] = Remaining$density[1]
            Peaks$level[p] = max(Remaining$level,na.rm=T)
            Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
            Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
            if(grepl("GC",sample)){
              if(Peaks$Start_var[p]>0.05){
                Peaks$Start_var[p] = 0.05
              }
              if(Peaks$End_var[p]>0.05){
                Peaks$End_var[p] = 0.05
              }
            } else {
              if(Peaks$Start_var[p]>0.1){
                Peaks$Start_var[p] = 0.1
              }
              if(Peaks$End_var[p]>0.1){
                Peaks$End_var[p] = 0.1
              }
            }
          } else {
            Peaks[p,] = Remaining[1,]
            Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
            Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
            Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
            Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
            Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
            Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
            Peaks$Height[p] = mean(Remaining$Height,na.rm = T)
            Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
            Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
            Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
            Peaks$density[p] = Remaining$density[1]
            Peaks$level[p] = max(Remaining$level,na.rm=T)
            Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
            Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
            if(grepl("GC",sample)){
              if(Peaks$Start_var[p]>0.05){
                Peaks$Start_var[p] = 0.05
              }
              if(Peaks$End_var[p]>0.05){
                Peaks$End_var[p] = 0.05
              }
            } else {
              if(Peaks$Start_var[p]>0.1){
                Peaks$Start_var[p] = 0.1
              }
              if(Peaks$End_var[p]>0.1){
                Peaks$End_var[p] = 0.1
              }
            }
          }
        } else {
          Peaks[p,] = Remaining[1,]
          Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
          Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
          Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
          Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
          Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
          Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
          Peaks$Height[p] = mean(Remaining$Height,na.rm = T)
          Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
          Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
          Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
          Peaks$density[p] = Remaining$density[1]
          Peaks$level[p] = max(Remaining$level,na.rm=T)
          Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
          Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
          if(grepl("GC",sample)){
            if(Peaks$Start_var[p]>0.05){
              Peaks$Start_var[p] = 0.05
            }
            if(Peaks$End_var[p]>0.05){
              Peaks$End_var[p] = 0.05
            }
          } else {
            if(Peaks$Start_var[p]>0.1){
              Peaks$Start_var[p] = 0.1
            }
            if(Peaks$End_var[p]>0.1){
              Peaks$End_var[p] = 0.1
            }
          }
        }
      }
      if(length(QualIons)==1){
        if(length(QualIon_peaks$Compound)==0){
          Peaks$Comment[p] = "qual not found"
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        } else {
          diff = 1:length(QuanIon_peaks$Compound)
          for(e in 1:length(QuanIon_peaks$Compound)){
            diff[e] = min(abs(QuanIon_peaks$RT[e]-QualIon_peaks$RT))
          }
          if(min(diff)<=0.1){
            Remaining = QuanIon_peaks[diff<=0.1,]
            aimed_time = compounds_without_qual$`retention time`[p]
            diffi = abs(Remaining$RT-aimed_time)
            Remaining = Remaining[diffi<=RT_range,]
            Remaining_RT_1 = plyr::round_any(Remaining$RT,0.05,f=floor)
            if(grepl("GC",sample)){
              Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 2)))
            } else {
              Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 1)))
            }
            
            allions_frame = as.data.frame(matrix(nrow=(nrow(QuanIon_peaks)+nrow(QualIon_peaks)),ncol=2))
            allions_frame$V1 = c(QuanIon_peaks$Compound,QualIon_peaks$Compound)
            if(grepl("GC",sample)){
              allions_frame$V2 = round(c(plyr::round_any(QuanIon_peaks$RT,0.05,f=floor),plyr::round_any(QualIon_peaks$RT,0.05,f=floor)),digits = 2)
            } else {
              allions_frame$V2 = round(c(plyr::round_any(QuanIon_peaks$RT,0.05,f=floor),plyr::round_any(QualIon_peaks$RT,0.05,f=floor)),digits = 1)
            }
            
            consensus_frame = as.data.frame(matrix(nrow=length(ions),ncol=(length(Remaining_RT$Var1)+1)))
            consensus_frame$V1 = ions
            colnames(consensus_frame) = c("ion",as.character(Remaining_RT$Var1))
            
            for(bu in 1:length(ions)){
              frameg = allions_frame[allions_frame$V1==ions[bu],]
              for(ub in 1:length(Remaining_RT$Var1)){
                if(grepl("GC",sample)){
                  consensus_frame[bu,(ub+1)] = (round(plyr::round_any(min(abs(as.numeric(as.character(Remaining_RT$Var1[ub]))-as.numeric(as.character(frameg$V2))),na.rm=T),0.05,f=floor),digits = 2)<=0.1)
                } else {
                  consensus_frame[bu,(ub+1)] = (round(plyr::round_any(min(abs(as.numeric(as.character(Remaining_RT$Var1[ub]))-as.numeric(as.character(frameg$V2))),na.rm=T),0.05,f=floor),digits = 1)<=0.1)
                }
              }
            }
            
            bingo_rows = c()
            
            for(bing in 1:length(Remaining_RT$Var1)){
              bingo_rows[bing] = length(consensus_frame[consensus_frame[,(bing+1)]==T,(bing+1)])
            }
            
            bingo_times = bingo_rows == length(ions)
            
            Remaining_RT = dplyr::filter(Remaining_RT,bingo_times==T)
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
            }
            
            if(nrow(Remaining_RT)==0){
              Peaks$Comment[p] = "no matched RTs"
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
            
            if(length(Remaining_RT$Var1)>1){
              Remaining_RT$intensities = 1:length(Remaining_RT$Var1)
              for(v in 1:length(Remaining_RT$Var1)){
                if(grepl("GC",sample)){
                  Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[v]],na.rm = T)
                } else {
                  Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[v]],na.rm = T)
                }
              }
              Remaining_RT = Remaining_RT[!is.nan(Remaining_RT$intensities),]
              if(grepl("GC",sample)){
                Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)-as.numeric(paste(Remaining_RT$Var1))),digits = 2)
              } else {
                Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)-as.numeric(paste(Remaining_RT$Var1))),digits = 1)
              }
              
              Remaining_RT$ratio_multiple_threshold = NA
              Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff<=inner_search_window_multiple_peaks] = 1/inner_intensity_ratio_multiple_peaks
              Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff>inner_search_window_multiple_peaks] = 1/outer_intensity_ratio_multiple_peaks
              
              Ratios = max(Remaining_RT$intensities)/Remaining_RT$intensities
              Rem_Ratios = Ratios[Ratios!=1]
              
              if(T%in%grepl("IS",ID)){
                Remaining_RT = Remaining_RT[Ratios<=(1/outer_intensity_ratio_multiple_peaks),]
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                }
              } else if(as.numeric(ID)%in%Multiple$ID){
                zu = match(as.numeric(ID),Multiple$ID)
                if(T%in%(Ratios>Remaining_RT$ratio_multiple_threshold)){
                  Remaining_RT = Remaining_RT[Ratios<=Remaining_RT$ratio_multiple_threshold,]
                  if(grepl("GC",sample)){
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
                  } else {
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
                  }
                }
                if(Multiple$multiple[zu]!=0){
                  if(grepl("GC",sample)){
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[Multiple$multiple[zu]],]
                  } else {
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[Multiple$multiple[zu]],]
                  }
                }
              } else if(mean(Rem_Ratios>=(1/outer_intensity_ratio_multiple_peaks),na.rm = T)==1){
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                }
              } else if((nrow(Remaining_RT[Remaining_RT$time_diff==min(Remaining_RT$time_diff),])>1)){
                Remaining_RT = Remaining_RT[Remaining_RT$time_diff<=min(Remaining_RT$time_diff,na.rm=T),]
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                }
              } else {
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                }
              }
            }
            
            if(length(Remaining_RT$Var1)==0){
              Peaks$Comment[p] = paste0(Peaks$Comment[p]," RTs do not match;")
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
            
            if(length(Remaining$Compound)==0){
              Peaks$Comment[p] = paste0(Peaks$Comment[p]," overlapping qual ions, but wrong RT;")
            }
            if(!is.na(zu)){
              if(Multiple$multiple[zu]==0){
                Peaks[p,] = Remaining[1,]
                Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
                Peaks$Start_RT[p] = min(Remaining$Start_RT,na.rm = T)
                Peaks$Start_RT_level[p] = min(Remaining$Start_RT_level,na.rm = T)
                Peaks$End_RT_level[p] = max(Remaining$End_RT_level,na.rm=T)
                Peaks$RT[p] = median(Remaining$RT,na.rm = T)
                Peaks$End_RT[p] = max(Remaining$End_RT,na.rm = T)
                Peaks$Height[p] = max(Remaining$Height,na.rm = T)
                Peaks$Width[p] = Peaks$End_RT[p]-Peaks$Start_RT[p]
                Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
                Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
                Peaks$density[p] = Remaining$density[1]
                Peaks$level[p] = max(Remaining$level,na.rm=T)
                Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
                Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
                if(grepl("GC",sample)){
                  if(Peaks$Start_var[p]>0.05){
                    Peaks$Start_var[p] = 0.05
                  }
                  if(Peaks$End_var[p]>0.05){
                    Peaks$End_var[p] = 0.05
                  }
                } else {
                  if(Peaks$Start_var[p]>0.1){
                    Peaks$Start_var[p] = 0.1
                  }
                  if(Peaks$End_var[p]>0.1){
                    Peaks$End_var[p] = 0.1
                  }
                }
              } else {
                Peaks[p,] = Remaining[1,]
                Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
                Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
                Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
                Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
                Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
                Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
                Peaks$Height[p] = max(Remaining$Height,na.rm = T)
                Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
                Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
                Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
                Peaks$density[p] = Remaining$density[1]
                Peaks$level[p] = max(Remaining$level,na.rm=T)
                Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
                Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
                if(grepl("GC",sample)){
                  if(Peaks$Start_var[p]>0.05){
                    Peaks$Start_var[p] = 0.05
                  }
                  if(Peaks$End_var[p]>0.05){
                    Peaks$End_var[p] = 0.05
                  }
                } else {
                  if(Peaks$Start_var[p]>0.1){
                    Peaks$Start_var[p] = 0.1
                  }
                  if(Peaks$End_var[p]>0.1){
                    Peaks$End_var[p] = 0.1
                  }
                }
              }
            } else {
              Peaks[p,] = Remaining[1,]
              Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
              Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
              Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
              Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
              Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
              Peaks$Height[p] = max(Remaining$Height,na.rm = T)
              Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
              Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
              Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
              Peaks$density[p] = Remaining$density[1]
              Peaks$level[p] = max(Remaining$level,na.rm=T)
              Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
              Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
              if(grepl("GC",sample)){
                if(Peaks$Start_var[p]>0.05){
                  Peaks$Start_var[p] = 0.05
                }
                if(Peaks$End_var[p]>0.05){
                  Peaks$End_var[p] = 0.05
                }
              } else {
                if(Peaks$Start_var[p]>0.1){
                  Peaks$Start_var[p] = 0.1
                }
                if(Peaks$End_var[p]>0.1){
                  Peaks$End_var[p] = 0.1
                }
              }
            }
          } else {
            Peaks$Comment[p] = "qual ion RT shift too large"
            final_ions = Peaks[p,]
            final_ions$ratio = 1
            final_ions$type = "quan"
            Peaks_list[[p]] = final_ions
          }
        }
      } 
      if(length(QualIons)>1){
        ionas = as.character(as.data.frame(table(QualIon_peaks$Compound))$Var1)
        if(length(ionas)==0){
          Peaks$Comment[p] = "qual ions not found"
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        }
        if(length(ionas)!=length(QualIons)&length(ionas)>0){
          #needed because if qual has NA, it would not be in QualIon_peaks anymore
          notfound = QualIons[!(QualIons%in%ionas)]
          Peaks$Comment[p] = paste0(paste0(notfound,collapse = ";")," not found")
          final_ions = Peaks[p,]
          final_ions$ratio = 1
          final_ions$type = "quan"
          Peaks_list[[p]] = final_ions
          next
        } else {
          mad = 1:length(QuanIon_peaks$Compound)
          for(qu in 1:length(QuanIon_peaks$Compound)){
            QualIon_peaks$diff = abs(QuanIon_peaks$RT[qu]-QualIon_peaks$RT)
            diff = 1:length(QualIons)
            for(qual in 1:length(QualIons)){
              diff[qual] = min(QualIon_peaks[QualIon_peaks$Compound==QualIons[qual],]$diff)
            }
            mad[qu] = mean(diff,na.rm = T)
          }
          if(min(mad)<=0.1){
            Remaining = QuanIon_peaks[mad<=0.1,]
            aimed_time = compounds_without_qual$`retention time`[p]
            diffi = abs(Remaining$RT-aimed_time)
            Remaining = Remaining[diffi<=RT_range,]
            Remaining_RT_1 = plyr::round_any(Remaining$RT,0.05,f=floor)
            if(grepl("GC",sample)){
              Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 2)))
            } else {
              Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 1)))
            }
            
            if(nrow(Remaining_RT)==0){
              Peaks$Comment[p] = "no matched RTs"
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
            
            allions_frame = as.data.frame(matrix(nrow=(nrow(QuanIon_peaks)+nrow(QualIon_peaks)),ncol=2))
            allions_frame$V1 = c(QuanIon_peaks$Compound,QualIon_peaks$Compound)
            if(grepl("GC",sample)){
              allions_frame$V2 = round(c(plyr::round_any(QuanIon_peaks$RT,0.05,f=floor),plyr::round_any(QualIon_peaks$RT,0.05,f=floor)),digits = 2)
            } else {
              allions_frame$V2 = round(c(plyr::round_any(QuanIon_peaks$RT,0.05,f=floor),plyr::round_any(QualIon_peaks$RT,0.05,f=floor)),digits = 1)
            }
            
            consensus_frame = as.data.frame(matrix(nrow=length(ions),ncol=(length(Remaining_RT$Var1)+1)))
            consensus_frame$V1 = ions
            colnames(consensus_frame) = c("ion",as.character(Remaining_RT$Var1))
            
            for(bu in 1:length(ions)){
              frameg = allions_frame[allions_frame$V1==ions[bu],]
              for(ub in 1:length(Remaining_RT$Var1)){
                if(grepl("GC",sample)){
                  consensus_frame[bu,(ub+1)] = (round(plyr::round_any(min(abs(as.numeric(as.character(Remaining_RT$Var1[ub]))-as.numeric(as.character(frameg$V2))),na.rm=T),0.05,f=floor),digits = 2)<=0.1)
                } else {
                  consensus_frame[bu,(ub+1)] = (round(plyr::round_any(min(abs(as.numeric(as.character(Remaining_RT$Var1[ub]))-as.numeric(as.character(frameg$V2))),na.rm=T),0.05,f=floor),digits = 1)<=0.1)
                }
              }
            }
            
            bingo_rows = c()
            
            for(bing in 1:length(Remaining_RT$Var1)){
              bingo_rows[bing] = length(consensus_frame[consensus_frame[,(bing+1)]==T,(bing+1)])
            }
            
            bingo_times = bingo_rows == length(ions)
            
            Remaining_RT = dplyr::filter(Remaining_RT,bingo_times==T)
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
            }
            
            if(length(Remaining_RT$Var1)>1){
              Remaining_RT$intensities = 1:length(Remaining_RT$Var1)
              for(v in 1:length(Remaining_RT$Var1)){
                if(grepl("GC",sample)){
                  Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[v]],na.rm = T)
                } else {
                  Remaining_RT$intensities[v] = max(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[v]],na.rm = T)
                }
              }
              Remaining_RT = Remaining_RT[!is.nan(Remaining_RT$intensities),]
              if(grepl("GC",sample)){
                Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)-as.numeric(paste(Remaining_RT$Var1))),digits = 2)
              } else {
                Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)-as.numeric(paste(Remaining_RT$Var1))),digits = 1)
              }
              
              Remaining_RT$ratio_multiple_threshold = NA
              Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff<=inner_search_window_multiple_peaks] = 1/inner_intensity_ratio_multiple_peaks
              Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff>inner_search_window_multiple_peaks] = 1/outer_intensity_ratio_multiple_peaks
              
              Ratios = max(Remaining_RT$intensities)/Remaining_RT$intensities
              Rem_Ratios = Ratios[Ratios!=1]
              
              if(as.numeric(ID)%in%Multiple$ID){
                zu = match(as.numeric(ID),Multiple$ID)
                if(T%in%(Ratios>Remaining_RT$ratio_multiple_threshold)){
                  Remaining_RT = Remaining_RT[Ratios<=Remaining_RT$ratio_multiple_threshold,]
                  if(grepl("GC",sample)){
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
                  } else {
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
                  }
                }
                if(Multiple$multiple[zu]!=0){
                  if(grepl("GC",sample)){
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[Multiple$multiple[zu]],]
                  } else {
                    Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[Multiple$multiple[zu]],]
                  }
                }
              } else if(mean(Rem_Ratios>=(1/outer_intensity_ratio_multiple_peaks),na.rm = T)==1){
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                }
              } else if((nrow(Remaining_RT[Remaining_RT$time_diff==min(Remaining_RT$time_diff,na.rm=T),])>1)){
                Remaining_RT = Remaining_RT[Remaining_RT$time_diff==min(Remaining_RT$time_diff,na.rm=T),]
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
                }
              } else {
                if(grepl("GC",sample)){
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                } else {
                  Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
                }
              }
            }
            
            if(length(Remaining_RT$Var1)==0){
              Peaks$Comment[p] = paste0(Peaks$Comment[p]," RTs do not match;")
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
            
            if(length(Remaining$Compound)==0){
              Peaks$Comment[p] = paste0(Peaks$Comment[p]," overlapping qual ions, but wrong RT;")
              final_ions = Peaks[p,]
              final_ions$ratio = 1
              final_ions$type = "quan"
              Peaks_list[[p]] = final_ions
              next
            }
            if(!is.na(zu)){
              if(Multiple$multiple[zu]==0){
                Peaks[p,] = Remaining[1,]
                Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
                Peaks$Start_RT[p] = min(Remaining$Start_RT,na.rm = T)
                Peaks$Start_RT_level[p] = min(Remaining$Start_RT_level,na.rm = T)
                Peaks$End_RT_level[p] = max(Remaining$End_RT_level,na.rm=T)
                Peaks$RT[p] = median(Remaining$RT,na.rm = T)
                Peaks$End_RT[p] = max(Remaining$End_RT,na.rm = T)
                Peaks$Height[p] = max(Remaining$Height,na.rm = T)
                Peaks$Width[p] = Peaks$End_RT[p]-Peaks$Start_RT[p]
                Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
                Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
                Peaks$density[p] = Remaining$density[1]
                Peaks$level[p] = max(Remaining$level,na.rm=T)
                Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
                Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
                if(grepl("GC",sample)){
                  if(Peaks$Start_var[p]>0.1){
                    Peaks$Start_var[p] = 0.1
                  }
                  if(Peaks$End_var[p]>0.1){
                    Peaks$End_var[p] = 0.1
                  }
                } else {
                  if(Peaks$Start_var[p]>0.1){
                    Peaks$Start_var[p] = 0.1
                  }
                  if(Peaks$End_var[p]>0.1){
                    Peaks$End_var[p] = 0.1
                  }
                }
              } else {
                Peaks[p,] = Remaining[1,]
                Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
                Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
                Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
                Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
                Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
                Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
                Peaks$Height[p] = max(Remaining$Height,na.rm = T)
                Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
                Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
                Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
                Peaks$density[p] = Remaining$density[1]
                Peaks$level[p] = max(Remaining$level,na.rm=T)
                Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
                Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
                if(grepl("GC",sample)){
                  if(Peaks$Start_var[p]>0.05){
                    Peaks$Start_var[p] = 0.05
                  }
                  if(Peaks$End_var[p]>0.05){
                    Peaks$End_var[p] = 0.05
                  }
                } else {
                  if(Peaks$Start_var[p]>0.1){
                    Peaks$Start_var[p] = 0.1
                  }
                  if(Peaks$End_var[p]>0.1){
                    Peaks$End_var[p] = 0.1
                  }
                }
              }
            } else {
              Peaks[p,] = Remaining[1,]
              Peaks$Sequence[p] = SAX_consensus(Remaining$Sequence)
              Peaks$Start_RT[p] = mean(Remaining$Start_RT,na.rm = T)
              Peaks$Start_RT_level[p] = mean(Remaining$Start_RT_level,na.rm=T)
              Peaks$End_RT_level[p] = mean(Remaining$End_RT_level,na.rm=T)
              Peaks$RT[p] = mean(Remaining$RT,na.rm = T)
              Peaks$End_RT[p] = mean(Remaining$End_RT,na.rm = T)
              Peaks$Height[p] = max(Remaining$Height,na.rm = T)
              Peaks$Width[p] = mean(Remaining$Width,na.rm = T)
              Peaks$Bg_Start[p] = min(Remaining$Bg_Start,na.rm = T)
              Peaks$Bg_End[p] = max(Remaining$Bg_End,na.rm = T)
              Peaks$density[p] = Remaining$density[1]
              Peaks$level[p] = max(Remaining$level,na.rm=T)
              Peaks$Start_var[p] = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
              Peaks$End_var[p] = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
              if(grepl("GC",sample)){
                if(Peaks$Start_var[p]>0.05){
                  Peaks$Start_var[p] = 0.05
                }
                if(Peaks$End_var[p]>0.05){
                  Peaks$End_var[p] = 0.05
                }
              } else {
                if(Peaks$Start_var[p]>0.1){
                  Peaks$Start_var[p] = 0.1
                }
                if(Peaks$End_var[p]>0.1){
                  Peaks$End_var[p] = 0.1
                }
              }
            }
          } else {
            Peaks$Comment[p] = "qual identified, but RT shift too high"
            final_ions = Peaks[p,]
            final_ions$ratio = 1
            final_ions$type = "quan"
            Peaks_list[[p]] = final_ions
            next
          }
        }
      }
    }
    if(length(QualIon_peaks$Compound[!is.na(QualIon_peaks$Height)])!=0){
      QualIon_peaks = QualIon_peaks[QualIon_peaks$RT>=(Peaks$RT[p]-0.1)&QualIon_peaks$RT<=(Peaks$RT[p]+0.1),]
      if(length(QualIon_peaks$Compound[!is.na(QualIon_peaks$Height)])!=0){
        confirming_ions = QualIon_peaks[1:length(as.data.frame(table(QualIon_peaks$Compound))$Var1),]
        ions = as.data.frame(table(QualIon_peaks$Compound))
        for(y in 1:length(confirming_ions$Compound)){
          confi = QualIon_peaks[QualIon_peaks$Compound==ions$Var1[y],]
          confirming_ions$Sequence[y] = SAX_consensus(confi$Sequence)
          confirming_ions$Compound[y] = confi$Compound[1]
          confirming_ions$mz[y] = confi$mz[1]
          confirming_ions$Comment[y] = confi$Comment[1]
          confirming_ions$Start_RT[y] = mean(confi$Start_RT,na.rm = T)
          confirming_ions$RT[y] = mean(confi$RT,na.rm = T)
          confirming_ions$End_RT[y] = mean(confi$End_RT,na.rm = T)
          confirming_ions$Sequence[y] = confi$Sequence[1]
          confirming_ions$Width[y] = mean(confi$Width,na.rm = T)
          confirming_ions$Height[y] = mean(confi$Height,na.rm = T)
          confirming_ions$Area[y] = mean(confi$Area,na.rm = T)
          confirming_ions$Start_var[y] = plyr::round_any((max(confirming_ions$Start_RT,na.rm = T)-min(confirming_ions$Start_RT,na.rm = T)),0.1,ceiling)
          confirming_ions$End_var[y] = plyr::round_any((max(confirming_ions$End_RT,na.rm = T)-min(confirming_ions$End_RT,na.rm = T)),0.1,ceiling)
        }
      } else {
        confirming_ions = NA
      }
    } else {
      confirming_ions = NA
    }
    final_ions = Peaks[p,]
    final_ions$ratio = 1
    final_ions$type = "quan"
    if(!is.na(confirming_ions)){
      confirming_ions$max_MINDIST = NA
      confirming_ions$ratio = confirming_ions$Height/Peaks$Height[p]
      if(T%in%(confirming_ions$ratio>1)){
        Peaks$Comment[p] = paste0(Peaks$Comment[p],"qual ion :",confirming_ions$mz[match(T,(confirming_ions$ratio>1))]," has higher intensity")
      }
      confirming_ions$type = "qual"
      confirming_ions$diff = NULL
      final_ions = rbind(final_ions,confirming_ions)
    }
    if(nrow(final_ions[!is.na(final_ions$RT),])>0){
      Consensus = SAX_consensus(final_ions$Sequence)
      for(vlu in 1:nrow(final_ions)){
        final_ions$max_MINDIST[vlu] = SAX_mindist(final_ions$Sequence[vlu],Consensus,final_ions$Nr_of_Points[vlu],SAX_reference_table)
      }
      Peaks$max_MINDIST[p] = max(final_ions$max_MINDIST,na.rm = T)
      Peaks$Nr_of_Points[p] = final_ions$Nr_of_Points[match(max(final_ions$max_MINDIST,na.rm = T),final_ions$max_MINDIST)]
      maximum_MINDIST = 0.5/sqrt(30)*sqrt(Peaks$Nr_of_Points[p])
      if(maximum_MINDIST>max_MINDIST){
        Peaks$Comment[p] = paste0(Peaks$Comment[p],"MINDIST > ",max_MINDIST,collapse = ";")
      }
      if(max(final_ions$max_MINDIST,na.rm = T)<maximum_MINDIST){
        Peaks$max_MINDIST[p] = maximum_MINDIST
      }
    }
    Peaks_list[[p]] = final_ions
    confirming_ions = NA
  }
  
  Peaks$Comment[is.na(Peaks$mz)] = paste0(Peaks$Comment[is.na(Peaks$mz)],"; Pattern not found")
  Peaks$mz[is.na(Peaks$mz)] = compounds_without_qual$`m/z`[is.na(Peaks$mz)]
  Peaks$Comment = stringr::str_remove_all(Peaks$Comment,"NA")
  
  ##############################################################################################
  #Check peaks on concentration/intensity-dependent level to see if overload shifts occur
  Intensity_dependent_shift_names = c()
  Intensity_dependent_shift = list()
  #Only checked for CAL since here a curve can be generated
  if(T%in%grepl("CAL",References_refined)){
    Concentration_levels = as.numeric(as.character(as.data.frame(table(Peaks_all$Concentration))$Var1))
    Concentration_levels = Concentration_levels[!is.na(Concentration_levels)]
    Concentration_levels = Concentration_levels[order(Concentration_levels,decreasing = F)]
  } else {
    Concentration_levels = References_refined
  }
  
  currentFunction = "check for intensity-dependent shift"
  
  if(T%in%grepl("CAL",References_refined)){
    for(p in 1:nrow(Peaks)){
      if(is.na(Peaks$RT[p])){
        next
      }
      aimed_time = Peaks$RT[p]
      intensity_dependent_frame = data.frame("int"=NA,"RT"=NA,"sequence"=NA,"concentration"=NA)
      empty = intensity_dependent_frame
      for(p_sub in 1:length(Concentration_levels)){
        zu = NA
        ID = substr(Peaks$Compound[p],1,4)
        selected = Peaks_all[grepl(ID,Peaks_all$Compound),]
        selected = selected[Concentration_levels[p_sub]==selected$Concentration,]
        ions = as.character(as.data.frame(table(selected$Compound))$Var1)
        
        QuanIon = ions[1]
        QuanIon_peaks = selected[selected$Compound%in%QuanIon,]
        QuanIon_peaks = QuanIon_peaks[!is.na(QuanIon_peaks$RT),]
        if(nrow(QuanIon_peaks)==0){
          intensity_dependent_frame = rbind(intensity_dependent_frame,empty)
          next
        }
        diffi = abs(QuanIon_peaks$RT-aimed_time)
        Remaining = QuanIon_peaks[diffi<=RT_range,]
        Remaining_RT_1 = plyr::round_any(Remaining$RT,0.05,f=floor)
        if(grepl("GC",sample)){
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 2)))
        } else {
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 1)))
        }
        
        if(length(Remaining_RT$Var)==0){
          intensity_dependent_frame = rbind(intensity_dependent_frame,empty)
          next
        }
        
        if(length(Remaining_RT$Var1)>1){
          Remaining_RT$intensities = 1:length(Remaining_RT$Var1)
          for(v in 1:length(Remaining_RT$Var1)){
            if(grepl("GC",sample)){
              Remaining_RT$intensities[v] = mean(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[v]],na.rm = T)
            } else {
              Remaining_RT$intensities[v] = mean(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[v]],na.rm = T)
            }
          }
          Remaining_RT = Remaining_RT[!is.nan(Remaining_RT$intensities),]
          if(grepl("GC",sample)){
            Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)-as.numeric(paste(Remaining_RT$Var1))),digits = 2)
          } else {
            Remaining_RT$time_diff = round(abs(round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)-as.numeric(paste(Remaining_RT$Var1))),digits = 1)
          }
          
          Remaining_RT$ratio_multiple_threshold = NA
          Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff<=inner_search_window_multiple_peaks] = 1/inner_intensity_ratio_multiple_peaks
          Remaining_RT$ratio_multiple_threshold[Remaining_RT$time_diff>inner_search_window_multiple_peaks] = 1/outer_intensity_ratio_multiple_peaks
          
          Ratios = max(Remaining_RT$intensities)/Remaining_RT$intensities
          Rem_Ratios = Ratios[Ratios!=1]
          
          if(as.numeric(ID)%in%Multiple$ID){
            zu = match(as.numeric(ID),Multiple$ID)
            if(T%in%(Ratios>Remaining_RT$ratio_multiple_threshold)){
              Remaining_RT = Remaining_RT[Ratios<=Remaining_RT$ratio_multiple_threshold,]
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)%in%Remaining_RT$Var1,]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)%in%Remaining_RT$Var1,]
              }
            }
            if(Multiple$multiple[zu]!=0){
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[Multiple$multiple[zu]],]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[Multiple$multiple[zu]],]
              }
            }
            #If 0, merge together
            Remaining[p,] = Remaining[1,]
            Remaining$Sequence = SAX_consensus(Remaining$Sequence)
            Remaining$Start_RT = min(Remaining$Start_RT,na.rm = T)
            Remaining$Start_RT_level = min(Remaining$Start_RT_level,na.rm = T)
            Remaining$End_RT_level = max(Remaining$End_RT_level,na.rm=T)
            Remaining$RT = mean(Remaining$RT,na.rm = T)
            Remaining$End_RT = max(Remaining$End_RT,na.rm = T)
            Remaining$Height = max(Remaining$Height,na.rm = T)
            Remaining$Width = Remaining$End_RT-Remaining$Start_RT
            Remaining$Bg_Start = min(Remaining$Bg_Start,na.rm = T)
            Remaining$Bg_End = max(Remaining$Bg_End,na.rm = T)
            Remaining$density = Remaining$density[1]
            Remaining$level = max(Remaining$level,na.rm=T)
            Remaining$Start_var = plyr::round_any((max(Remaining$Start_RT_level,na.rm = T)-min(Remaining$Start_RT_level,na.rm = T)),0.1,ceiling)
            Remaining$End_var = plyr::round_any((max(Remaining$End_RT_level,na.rm = T)-min(Remaining$End_RT_level,na.rm = T)),0.1,ceiling)
          } else if(mean(Rem_Ratios>=(1/outer_intensity_ratio_multiple_peaks),na.rm = T)==1){
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            } else {
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            }
          } else if((nrow(Remaining_RT[Remaining_RT$time_diff>=min(Remaining_RT$time_diff),]>1))){
            Remaining_RT = Remaining_RT[Remaining_RT$time_diff<=min(Remaining_RT$time_diff,na.rm=T),]
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
            } else {
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(max(Remaining_RT$intensities),Remaining_RT$intensities)],]
              }
            }
          } else {
            if(grepl("GC",sample)){
              Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
            } else {
              if(grepl("GC",sample)){
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
              } else {
                Remaining = Remaining[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[match(min(Remaining_RT$time_diff),Remaining_RT$time_diff)],]
              }
            }
          }
        }
        intensity_dependent_frame = rbind(intensity_dependent_frame,empty)
        intensity_dependent_frame$int[p_sub] = Remaining$Height
        intensity_dependent_frame$RT[p_sub] = Remaining$RT
        intensity_dependent_frame$sequence[p_sub] = Remaining$Sequence
        intensity_dependent_frame$concentration[p_sub] = Concentration_levels[p_sub]
      }
      intensity_dependent_frame = intensity_dependent_frame[!is.na(intensity_dependent_frame$RT),]
      
      if(0.0%in%intensity_dependent_frame$concentration){
        p_sub = match(0.0,Concentration_levels)
        zu = NA
        ID = substr(Peaks$Compound[p],1,4)
        selected = Peaks_all[grepl(ID,Peaks_all$Compound),]
        selected = selected[Concentration_levels[p_sub]==selected$Concentration,]
        ions = as.character(as.data.frame(table(selected$Compound))$Var1)
        
        QuanIon = ions[1]
        QuanIon_peaks = selected[selected$Compound%in%QuanIon,]
        QuanIon_peaks = QuanIon_peaks[!is.na(QuanIon_peaks$RT),]
        if(nrow(QuanIon_peaks)==0){
          intensity_dependent_frame = rbind(intensity_dependent_frame,empty)
          next
        }
        diffi = abs(QuanIon_peaks$RT-aimed_time)
        Remaining = QuanIon_peaks[diffi<=RT_range,]
        Remaining_RT_1 = plyr::round_any(Remaining$RT,0.05,f=floor)
        if(grepl("GC",sample)){
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 2)))
        } else {
          Remaining_RT = as.data.frame(table(round(Remaining_RT_1,digits = 1)))
        }
        
        if(length(Remaining_RT$Var1)>1){
          Remaining_RT$intensities = 1:length(Remaining_RT$Var1)
          for(v in 1:length(Remaining_RT$Var1)){
            if(grepl("GC",sample)){
              Remaining_RT$intensities[v] = mean(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 2)==Remaining_RT$Var1[v]],na.rm = T)
            } else {
              Remaining_RT$intensities[v] = mean(Remaining$Height[round(plyr::round_any(Remaining$RT,0.05,f=floor),digits = 1)==Remaining_RT$Var1[v]],na.rm = T)
            }
          }
          Potential_interfering = Remaining_RT[!is.nan(Remaining_RT$intensities),]
          colnames(Potential_interfering) = c("RT","sequence","int")
          Potential_interfering$concentration = 0
          Potential_interfering$RT = as.numeric(as.character(Potential_interfering$RT))
        } else {
          Potential_interfering = empty
        }
      } else {
        Potential_interfering = empty
      }
      
      if(nrow(Potential_interfering)>1){
        for(inf in 1:nrow(Potential_interfering)){
          if(grepl("GC",sample)){
            diffs = round(abs(Potential_interfering$RT[inf]-round(plyr::round_any(intensity_dependent_frame$RT,0.05,floor),digits=2)),digits = 1)
          } else {
            diffs = round(abs(Potential_interfering$RT[inf]-round(plyr::round_any(intensity_dependent_frame$RT,0.05,floor),digits=1)),digits = 1)
          }
          Potential_interfering_frame = intensity_dependent_frame[diffs<=0.1,]
          model = tryCatch(lm(int~concentration,data=Potential_interfering_frame),
                           error=function(e){model=NA})
          if(is.na(model)|nrow(Potential_interfering_frame)<4){
            intensity_dependent_frame = intensity_dependent_frame[diffs>0.1,]
            R2 = 0
          } else {
            R2 = summary(model)$r.squared
          }
          if(R2 < 0.7){
            intensity_dependent_frame = intensity_dependent_frame[diffs>0.1,]
          }
        }
      } 
      intensity_dependent_frame = intensity_dependent_frame[!is.na(intensity_dependent_frame$RT),]
      
      if(nrow(intensity_dependent_frame)<4){
        next
      }
      if(grepl("GC",sample)){
        aimed_time1 = round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)
        potential_times = as.numeric(as.character(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 2)))
      } else {
        aimed_time1 = round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)
        potential_times = as.numeric(as.character(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 1)))
      }
      diff_times = abs(potential_times-aimed_time1)
      intensity_dependent_frame = intensity_dependent_frame[diff_times<=0.5,]
      diff_times = diff_times[diff_times<=0.5]
      test = tryCatch(outliers::grubbs.test(diff_times,type=10,two.sided = T),
                      error=function(e){test=data.frame("p.value"=Inf)})
      if(test$p.value<0.05){
        res <- str_match(test$alternative, "value\\s*(.*?)\\s*is")
        Outlier = as.numeric(res[,2])
        intensity_dependent_frame = intensity_dependent_frame[intensity_dependent_frame$RT!=Outlier,]
      }
      model = tryCatch(lm(int~concentration,data=intensity_dependent_frame),
                       error=function(e){model=NA})
      if(is.na(model)){
        next
      }
      residuals = summary(model)$residuals
      residual_error = mean(abs(residuals/intensity_dependent_frame$int*100))
      R2 = summary(model)$r.squared
      model2 = tryCatch(lm(RT~concentration,data=intensity_dependent_frame),
                        error=function(e){model2=NA})
      if(is.na(model2)){
        next
      }
      R2_2 = summary(model2)$r.squared
      if(grepl("GC",sample)){
        aimed_time1 = round(plyr::round_any(aimed_time,0.05,f=floor),digits = 2)
        potential_times = as.numeric(as.character(as.data.frame(table(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 2)))$Var1))
        potential_frequencies = as.numeric(as.character(as.data.frame(table(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 2)))$Freq))
      } else {
        aimed_time1 = round(plyr::round_any(aimed_time,0.05,f=floor),digits = 1)
        potential_times = as.numeric(as.character(as.data.frame(table(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 1)))$Var1))
        potential_frequencies = as.numeric(as.character(as.data.frame(table(round(plyr::round_any(intensity_dependent_frame$RT,0.05,f=floor),digits = 1)))$Freq))
      }
      diff_times = abs(potential_times-aimed_time1)
      if(grepl("GC",sample)){
        sum_frequencies = sum(potential_frequencies[diff_times<=0.05],na.rm = T)
      } else {
        sum_frequencies = sum(potential_frequencies[diff_times<=0.1],na.rm = T)
      }
      if(length(potential_times)==1){
        ratio_same_RT = 1
      } else {
        ratio_same_RT = sum_frequencies/nrow(intensity_dependent_frame)
      }
      
      if(R2 > 0.7&R2_2>0.7&ratio_same_RT<=0.5&residual_error<=20){
        intensity_dependent_shift = data.frame("int"=rep(NA,nrow(intensity_dependent_frame)),
                                               "shift"=rep(NA,nrow(intensity_dependent_frame)))
        intensity_dependent_shift$int = plyr::round_any(as.numeric(intensity_dependent_frame$int),accuracy = 1)
        intensity_dependent_shift$shift = intensity_dependent_frame$RT-aimed_time
        intensity_dependent_shift = intensity_dependent_shift[order(intensity_dependent_shift$int,decreasing = F),]
        Peaks$Sequence[p] = SAX_consensus(intensity_dependent_frame$sequence)
        mindists = rep(NA,nrow(intensity_dependent_frame))
        for(min in 1:nrow(intensity_dependent_frame)){
          mindists[min] = SAX_mindist(Peaks$Sequence[p],intensity_dependent_frame$sequence[min],
                                      Peaks$Nr_of_Points[p],SAX_reference_table)
        }
        Peaks$max_MINDIST[p] = max(mindists,na.rm = T)
        Intensity_dependent_shift_names = append(Intensity_dependent_shift_names,Remove_ion(Peaks$Compound[p]))
        Intensity_dependent_shift = c(Intensity_dependent_shift,list(intensity_dependent_shift))
      } else if(R2 > 0.7){
        Peaks$Sequence[p] = SAX_consensus(intensity_dependent_frame$sequence)
        mindists = rep(NA,nrow(intensity_dependent_frame))
        for(min in 1:nrow(intensity_dependent_frame)){
          mindists[min] = SAX_mindist(Peaks$Sequence[p],intensity_dependent_frame$sequence[min],
                                      Peaks$Nr_of_Points[p],SAX_reference_table)
        }
        Peaks$max_MINDIST[p] = max(mindists,na.rm = T)
      } else if(T%in%grepl("IS",ID)|T%in%grepl("QA",References_refined)){
        Peaks$max_MINDIST[p] = Inf #internal standards have very unique masses and no calibration -> Ignore Mindist
      }
    }
  } else {
    for(p in 1:nrow(Peaks)){
      if(is.na(Peaks$RT[p])){
        next
      }
      aimed_time = Peaks$RT[p]
      intensity_dependent_frame = data.frame("int"=NA,"RT"=NA,"sequence"=NA)
      empty = intensity_dependent_frame
      
    }
  }
  
  names(Intensity_dependent_shift) = Intensity_dependent_shift_names
  Intensity_dependent_shift = Intensity_dependent_shift[!is.na(Intensity_dependent_shift_names)]
  
  Peaks$Compound = sapply(compounds_without_qual$identity,Remove_ion)
  Peaks$Start_var = as.numeric(Peaks$Start_var)
  Peaks$End_var = as.numeric(Peaks$End_var)
  
  Peaks$Start_var =  as.numeric(plyr::round_any(Peaks$Width,f=ceiling,accuracy=0.05)*1.5)
  Peaks$End_var =  as.numeric(plyr::round_any(Peaks$Width,f=ceiling,accuracy=0.05)*1.5)
  Peaks$RT_tol =  as.numeric(plyr::round_any(Peaks$Width,f=ceiling,accuracy=0.05))
  Peaks$Start_var[is.na(Peaks$Start_var)] = maximum_RT_tolerance_start
  Peaks$End_var[is.na(Peaks$End_var)] = maximum_RT_tolerance_end
  Peaks$RT_tol[is.na(Peaks$RT_tol)] = maximum_RT_tolerance
  
  Peaks$Start_var[Peaks$Start_var>maximum_RT_tolerance_start] = maximum_RT_tolerance_start
  Peaks$End_var[Peaks$End_var>maximum_RT_tolerance_end] = maximum_RT_tolerance_end
  Peaks$RT_tol[Peaks$RT_tol>maximum_RT_tolerance] = maximum_RT_tolerance
  
  Peaks$max_MINDIST[Peaks$max_MINDIST<max_MINDIST] = max_MINDIST
  
  Peaks_list_IS = Peaks_list[sapply(Peaks$Compound,extract_first_two)=="IS"]
  Peaks_IS = Peaks[sapply(Peaks$Compound,extract_first_two)=="IS",]
  
  Output = list()
  Output$Peaks = Peaks
  Output$Peaks_all = Peaks_all
  Output$References_refined = References_refined
  Output$Peaks_list = Peaks_list
  Output$Peaks_IS = Peaks_IS
  Output$Peaks_list_IS = Peaks_list_IS
  Output$Intensity_dependent_shift = Intensity_dependent_shift
  Output$compounds_shift_corrected = compounds_shift_corrected
  Output$compounds_without_qual = compounds_without_qual
  return(Output)
}

#'calculate_IS_dependent_shift
#'
#'Checks internal standards in each sample and calculates potential shifts
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @export
calculate_IS_dependent_shift = function(sample,set.env=F){
  if(set.env==T){
    environment(generate_peaklist_and_calculate_intensity_dependent_shift) = .GlobalEnv
    list2env(generate_peaklist_and_calculate_intensity_dependent_shift(sample,set.env=T),envir = .GlobalEnv)
  }
  
  n_iter <- length(Peaks_IS$Compound)
  
  pb <- txtProgressBar(min = 0,
                       max = n_iter,
                       style = 3,
                       width = n_iter, 
                       char = "=") 
  
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  currentFunction = "IS_analysis"
  
  Peaklists_IS = list()
  length(Peaklists_IS) = length(files)
  
  cl <- makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  environment(errors) = environment()
  tryCatch({
    for(j in 1:length(files)){
      print(paste0("Analysis of internal standards: ",files[j]))
      
      Peaklist_file_IS = foreach(t=1:length(Peaks_IS$Compound),
                                 .packages = c("seewave","dplyr","rawrr","RaMS","zoo","stats","MESS","tcpl","stringr"),
                                 .combine = "rbind",
                                 .multicombine = T,
                                 .options.snow=opts,
                                 .export = c("sample_path","files","sample","Multiple","RT_range","SAX_reference_table","compounds_shift_corrected","Solvent_Blank_ID","Peaks_IS","Peaks_list_IS","FullscanMS1","method_time","ppm_val","functions",names(functions),"j","alphabet_size","zigzag_trigger_threshold","extended_baseline_factor","normal_background_quantile","width_factor_background","sample_search_window","higher_background_quantile","minimum_background_ratio","width_smoothing","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","maximum_peak_width","minimum_peak_area","minimum_nr_of_datapoints_per_peak")) %dopar% {
                                   IS_analysis(sample_path,sample,SAX_reference_table,Solvent_Blank_ID,Multiple,RT_range,files,Peaks=Peaks_IS,Peaks_list=Peaks_list_IS,compounds_shift_corrected,FullscanMS1,Scan_filters,ppm_val,method_time,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,minimum_background_ratio,extended_baseline_factor,width_smoothing,width_factor_background,sample_search_window,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,m,t,j)
                                 }
      
      Peaklists_IS[[j]] = Peaklist_file_IS
      
    }
  },error = errors)
  
  close(pb)
  stopCluster(cl)
  
  IS_dependent_shift = as.data.frame(matrix(nrow=nrow(Peaks_IS),ncol=length(files)+1))
  colnames(IS_dependent_shift) = c("ISTD",files)
  IS_dependent_shift$ISTD = Peaks_IS$Compound
  
  expected_times = data.frame("ISTD"=Peaks_IS$Compound,"RT"=Peaks_IS$RT)
  expected_times$RT = as.numeric(expected_times$RT)
  for(r in 1:length(files)){
    for(si in 1:length(Peaks_IS$Compound)){
      asi = (as.numeric(Peaklists_IS[[r]]$RT[si])-expected_times$RT[si])
      if(is.na(asi)){
        IS_dependent_shift[si,(r+1)] = 0
      } else {
        IS_dependent_shift[si,(r+1)] = asi
      }
    }
  }
  
  Output = list()
  Output$IS_dependent_shift = IS_dependent_shift
  return(Output)
}

#'target_screening_of_files
#'
#'Performs the target screening for each file within the sample batch
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @param gen.plots boolean if plots should be generated
#' @param use.MINDIST boolean if MINDIST should be applied as peak selection criterion
#' @param use.area boolean, if TRUE uses peak area, otherwise peak height
#' @export
target_screening_of_files = function(sample,set.env=F,gen.plots,use.MINDIST,use.area){
  if(set.env==T){
    environment(calculate_IS_dependent_shift) = .GlobalEnv
    list2env(calculate_IS_dependent_shift(sample,set.env=T),envir = .GlobalEnv)
  }
  
  Results = as.data.frame(matrix(nrow = length(Peaks$Compound),
                                 ncol = length(files)+15))
  
  dir.create(results_path)
  if(T%in%grepl(".mzML",References)){
    colnames(Results) = c("Compound","Comment","mz","RT_start","RT","RT_end","ISTD","CV IS","Valid references","Valid calibration","R2_calibration","Min_Calibration","Max_Calibration","Matrix correction","Blank values",str_remove(files,".mzML"))
  } else {
    colnames(Results) = c("Compound","Comment","mz","RT_start","RT","RT_end","ISTD","CV IS","Valid references","Valid calibration","R2_calibration","Min_Calibration","Max_Calibration","Matrix correction","Blank values",str_remove(files,".raw"))
  }
  
  n_iter <- length(files)
  
  pb <- txtProgressBar(min = 0,
                       max = n_iter,
                       style = 3,
                       width = n_iter, 
                       char = "=") 
  
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  Peaklists = list()
  length(Peaklists) = length(files)
  cl <- makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  
  currentFunction = "sample_analysis"
  environment(errors) = environment()
  tryCatch({
    for(j in 1:length(files)){
      init[j] <- Sys.time()
      
      Peaklist_file = foreach(t=1:length(Peaks$Compound),
                              .packages = c("seewave","dplyr","rawrr","zoo","stats","MESS","tcpl","stringr","ggplot2","ggtext","RaMS"),
                              .combine = "rbind",
                              .multicombine = T,
                              .export = c("sample_path","results_path","sample","files","Peaks","Peaks_list","Solvent_Blank_ID","IS_dependent_shift","Intensity_dependent_shift","IS_Assignment","compounds_shift_corrected","RT_range","FullscanMS1","precursor_window","inner_intensity_ratio_multiple_peaks","Scan_filters","Multiple","ppm_val","method_time","gen.plots","use.MINDIST","SAX_reference_table","IS_deviation","functions",names(functions),"j","alphabet_size","zigzag_trigger_threshold","normal_background_quantile","higher_background_quantile","minimum_background_ratio","extended_baseline_factor","width_smoothing","width_factor_background","sample_search_window","maximum_allowed_shift","maximum_allowed_shift_ratio","maximum_nr_of_a","minimum_nr_of_high_intensity_letters","minimum_nr_of_datapoints_per_peak","minimum_peak_area","maximum_peak_width","minimum_cutoff_intensity_factor","minimum_confirming_peak_height","max_MINDIST","minimum_datapoints_per_sample_peak","minimum_qualitative_threshold")) %dopar% {
                                sample_analysis(sample_path,results_path,sample,files,Peaks,Peaks_list,Solvent_Blank_ID,IS_dependent_shift,Intensity_dependent_shift,Multiple,FullscanMS1,precursor_window,inner_intensity_ratio_multiple_peaks,Scan_filters,IS_Assignment,compounds_shift_corrected,RT_range,ppm_val,method_time,IS_deviation,gen.plots,use.MINDIST,SAX_reference_table,alphabet_size,zigzag_trigger_threshold,normal_background_quantile,higher_background_quantile,minimum_background_ratio,extended_baseline_factor,width_smoothing,width_factor_background,sample_search_window,maximum_nr_of_a,minimum_nr_of_high_intensity_letters,maximum_peak_width,minimum_peak_area,minimum_nr_of_datapoints_per_peak,maximum_allowed_shift,maximum_allowed_shift_ratio,minimum_cutoff_intensity_factor,minimum_confirming_peak_height,max_MINDIST,minimum_datapoints_per_sample_peak,minimum_qualitative_threshold,t,j,rep=F)
                              }
      
      Peaklists[[j]] = Peaklist_file
      
      end[j] <- Sys.time()
      
      setTxtProgressBar(pb, j)
      time <- round(seconds_to_period(sum(end - init)), 0)
      
      est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
      remainining <- round(seconds_to_period(est), 0)
      
      cat(paste(" Finished: ",files[j],"\n",
                " // Execution time:", time,
                " // Estimated time remaining:", remainining), "")
    }
  },error=errors)
  
  stopCluster(cl)
  
  Results$Compound = Peaks$Compound
  Results$mz = Peaks$mz
  Results$Comment = Peaks$Comment
  Results$RT_start = Peaks$Start_RT
  Results$RT = Peaks$RT
  Results$RT_end = Peaks$End_RT
  
  if(use.area==T){
    for(h in 1:length(Peaklists)){
      Results[,h+15] = Peaklists[[h]]$Area
    }
  } else {
    for(h in 1:length(Peaklists)){
      Results[,h+15] = Peaklists[[h]]$Height
    }
  }
  
  Results_RT = Results
  
  for(h in 1:length(Peaklists)){
    Results_RT[,h+15] = Peaklists[[h]]$RT
  }
  
  currentFunction = "Summarize_IS"
  
  names(Peaklists) = files
  files_wo_Solvent_Blank = files[!grepl(Solvent_Blank_ID,files)]
  Peaklists_wo_Solvent_Blank = Peaklists[files%in%files_wo_Solvent_Blank]
  IS_l = Peaks[Peaks$mz%in%IS$`m/z`,]
  IS_values = as.data.frame(matrix(nrow = length(IS_l$Compound),ncol = length(files_wo_Solvent_Blank)+1))
  if(T%in%grepl(".mzML",References)){
    colnames(IS_values) = c("IS",str_remove(files_wo_Solvent_Blank,".mzML"))
  } else {
    colnames(IS_values) = c("IS",str_remove(files_wo_Solvent_Blank,".raw"))
  }
  list = Peaklists_wo_Solvent_Blank[[1]]
  iss = list[list$mz%in%IS$`m/z`,]
  IS_values$IS = iss$Compound
  if(use.area==T){
    IS_values[,2] = iss$Area
  } else {
    IS_values[,2] = iss$Height
  }
  
  if(length(Peaklists_wo_Solvent_Blank)>1){
    for(is in 2:length(Peaklists_wo_Solvent_Blank)){
      list = Peaklists_wo_Solvent_Blank[[is]]
      iss = list[list$mz%in%IS$`m/z`,]
      if(use.area==T){
        IS_values[,is+1] = iss$Area
      } else {
        IS_values[,is+1] = iss$Height
      }
      
    }
  }
  
  IS_values$Mean = rowMeans(IS_values[,2:length(IS_values)],na.rm = T)
  IS_values$CV = NA
  for(w in 1:length(IS_values$IS)){
    IS_values$CV[w] = (sd(IS_values[w,2:length(IS_values)-2],na.rm = T)/IS_values$Mean[w])*100
  }
  
  currentFunction = "IS_metric_plots"
  for(i in 1:length(Results$Compound)){
    ID = substr(Results$Compound[i],1,4)
    if(T %in% grepl(ID,IS$ID)){
      ISTD_ID = ID
    } else {
      ISTD = IS_Assignment[grepl(ID,IS_Assignment$Compound),]$ISTD
      ISTD_ID = substr(ISTD,1,4)
    }
    u = match(Results[grepl(ISTD_ID,Results$Compound),]$Compound,Results$Compound)
    Results$ISTD[i] = Results$Compound[u]
    is = Results$ISTD[i]
    hm = grepl(ISTD_ID,IS_values$IS)
    k = match(T,hm)
    Results$`CV IS`[i] = IS_values$CV[k]
  }
  
  PDF_path = paste0(results_path,"/",Sys.Date(),"_IS_metric_plots.pdf")
  pdf(file=PDF_path) 
  
  for(is in 1:length(IS_values$IS)){
    levels = as.numeric(IS_values[is,2:(length(IS_values)-2)])
    if(mean(is.na(levels))==1){
      next
    }
    deviation = IS_deviation*mean(levels)
    
    # Use CV as measure of variability
    cv <- 100 * (sd(levels,na.rm=T) / mean(levels,na.rm = T))
    
    plot_title <- paste(IS_values$IS[is], "-", "CV:", round(cv, 2), "%")
    
    if(use.area==T){
      label = "Peak Area"
    } else {
      label = "Peak Height"
    }
    
    # Calculate the y-axis limits
    y_range <- range(levels)
    y_diff <- diff(y_range)
    deviation_limit <- ifelse(deviation > y_diff, 2 * deviation, 2 * y_diff)
    y_min <- min(levels) - deviation_limit
    y_max <- max(levels) + deviation_limit
    
    # Generate the plot
    p <- ggplot(data = NULL, aes(x = factor(seq_along(colnames(IS_values)[2:(length(IS_values)-2)])), y = levels)) +
      geom_point(shape = 16, size = 4, color = "black") +
      geom_hline(yintercept = mean(levels,na.rm=T), linetype = "dashed", color = "red") +
      geom_hline(yintercept = mean(levels,na.rm = T) + deviation, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = mean(levels,na.rm = T) - deviation, linetype = "dashed", color = "orange") +
      #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = mean(levels,na.rm = T) - deviation, ymax = mean(levels,na.rm = T) + deviation),
      #          fill = "lightgreen", alpha = 0.05) +
      labs(x = NULL, y = label, title = plot_title) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
            axis.text.y = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      scale_x_discrete(labels = colnames(IS_values)[2:(length(IS_values)-2)]) +
      coord_cartesian(ylim = c(y_min, y_max))
    
    print(p)
    
  }
  dev.off()
  
  Output = list()
  Output$Results = Results
  Output$Results_RT = Results_RT
  Output$Peaklists = Peaklists
  Output$IS_values = IS_values
  
  return(Output)
  
}

#'quantify_samples
#'
#'Performs quantitative analysis (recovery/quantitation) of all samples
#'
#' @param sample the name of the sample / file batch
#' @param set.env boolean if the environment of the R session should be setup from scratch
#' @param max_calibration to what extent (ratio, i.e. 1.2) are values above the maximum concentration still considered valid
#' @param use.area boolean, if TRUE uses peak area, otherwise peak height
#' @param quantified_unit the unit of the samples and final unit of quantified values
#' @param IS_deviation the value of accepted deviation for the IS in fractions (i.e. 0.2)
#' @param use.MINDIST boolean if MINDIST should be applied as peak selection criterion
#' @param gen.plots boolean if plots should be generated
#' @param Reference_deviation the value of accepted deviation for the IS in fractions (i.e. 0.2)
#' @param minimum_r2 the lowest accepted coefficient of determination for valid calibrations
#' @param maximum_error_of_quantitation the maximal accepted deviation in % of known and calculated concentrations
#' @param minimum_nr_of_datapoints_calibration the lowest number of data points a valid calibration should consist of
#' @export
quantify_samples = function(sample,set.env=F,max_calibration,use.area,quantified_unit,IS_deviation,use.MINDIST,gen.plots,Reference_deviation,minimum_r2,maximum_error_of_quantitation,minimum_nr_of_datapoints_calibration){
  
  if(set.env==T){
    environment(target_screening_of_files) = .GlobalEnv
    list2env(target_screening_of_files(sample,set.env=T,gen.plots = T,use.MINDIST = use.MINDIST,use.area=use.area,IS_deviation = IS_deviation))
  }
  
  currentFunction = "blank_effects"
  
  #Check if there are Peaks frequently found in solvent Blanks (without IS of course) - if yes, subtract from all other signals since system-dependent (not sample-dependent)
  Blanks = Results %>% dplyr::select(contains(c("Compound",Solvent_Blank_ID)))
  for(i in 1:length(Blanks$Compound)){
    Blanks[i,c(F,is.na(Blanks[i,2:length(Blanks)]))] = 0
    Too_high = Blanks[i,2:length(Blanks)] > 0
    if(length(Too_high[Too_high==T])>=0.33*(length(Blanks)-1)){
      Results$`Blank values`[i] = median(as.numeric(Blanks[i,c(F,Too_high)]),na.rm = T)
      Results[i,16:length(Results)] = Results[i,16:length(Results)]-Results$`Blank values`[i]
      Results[i,c(rep(F,15),is.na(Results[i,16:length(Results)]))] = 0
      Results[i,c(rep(F,15),Results[i,16:length(Results)]<=0)] = NA
    } else {
      Results$`Blank values`[i] = 0
    }
  }
  
  #Get a qualitative summary for all compounds in all samples
  currentFunction = "results_qualitative"
  
  Results_qualitative = Results[,c(1,16:length(Results))]
  
  Results_qualitative[1:length(Results_qualitative$Compound),2:length(Results_qualitative)] = NA
  
  for(w in 1:length(Peaks$Compound)){
    for(wj in 1:length(Peaklists)){
      if(!is.na(Peaklists[[wj]]$Height[w])){
        if(Peaklists[[wj]]$Height[w]>Peaklists[[wj]]$LOD[w]&(Peaklists[[wj]]$Area[w]-Results$`Blank values`[w])>0){
          Results_qualitative[w,wj+1] = T
        } else {
          Results_qualitative[w,wj+1] = F
        }
      } else {
        Results_qualitative[w,wj+1] = F
      }
    }
    
    VAlues = Results[w,16:length(Results)]
    
    VAlues[is.na(VAlues)] = Inf
    
    Results_qualitative[w,c(F,VAlues<(Results$`Blank values`[w]*2))] = "Background" #x2 because 1x subtracted and I want at least 3x to be sure that valid peak
    
  }
  
  for(w in 1:length(Peaks$Compound)){
    Results_qualitative[w,as.character(Results_qualitative[w,])=="0"] = "FALSE"
    Results_qualitative[w,as.character(Results_qualitative[w,])=="1"] = "TRUE"
  }
  
  #Only accept Peaks which are above LOQ - only if quantified
  if((T%in%grepl("CAL",References_refined))){
    for(w in 1:length(Peaks$Compound)){
      for(wj in 1:length(Peaklists)){
        LOQs = Peaklists[[wj]]$LOQ[w]
        Heights = Peaklists[[wj]]$Height[w]
        if(is.na(Heights)){
          next
        } else {
          if(Heights<LOQs){
            Results[w,wj+15] = NA
          }
        }
      }
    }
  }
  
  #Compare Values of solvent-matched reference to matrix-matched reference
  #Adjust ratios according to matrix effects, blow-down or blank effects and add IS info
  if(T%in%grepl(".mzML",files)){
    Reference_Values = Results %>% dplyr::select(contains(c("Compound","ISTD",str_remove(References,".mzML"))))
  } else {
    Reference_Values = Results %>% dplyr::select(contains(c("Compound","ISTD",str_remove(References,".raw"))))
  }
  currentFunction = "Matrix"
  Matrix_Values = Results %>% dplyr::select(contains(c("Compound","matrix","Matrix")))
  
  if(length(Matrix_Values)>2){
    for(i in 1:length(Matrix_Values$Compound)){
      k = match(Results$ISTD[i],Results$Compound)
      if(ncol(Matrix_Values)>3){
        matrix = rowMeans(Matrix_Values[i,3:length(Matrix_Values)],na.rm = T)
        matrix_IS = rowMeans(Matrix_Values[k,3:length(Matrix_Values)],na.rm = T)
        reference = rowMeans(Reference_Values[i,3:length(Reference_Values)],na.rm = T)
        reference_IS = rowMeans(Reference_Values[k,3:length(Reference_Values)],na.rm = T)
      } else {
        matrix = Matrix_Values[i,3]
        matrix_IS = Matrix_Values[k,3]
        reference = rowMeans(Reference_Values[i,3:length(Reference_Values)],na.rm = T)
        reference_IS = rowMeans(Reference_Values[k,3:length(Reference_Values)],na.rm = T)
      }
      matrix_effect = matrix/reference
      matrix_IS_effect = matrix_IS/reference_IS
      matrix_effect = matrix_effect*matrix_IS_effect
      if(!is.nan(matrix_effect)&!is.na(matrix_effect)){
        patterns = paste0("matrix","Matrix",Reference_ID,Solvent_Blank_ID,collapse = "|")
        not_matrix = !(grepl(patterns,colnames(Results)))
        not_matrix = c(rep(F,15),not_matrix[16:length(not_matrix)])
        Results[i,not_matrix] = Results[i,not_matrix]/matrix_effect
        Results$`Matrix correction`[i] = matrix_effect
      } else {
        Results$`Matrix correction`[i] = 0
      }
    }
  }
  
  currentFunction = "Results_ratio"
  Results_ratio = Results
  
  for(r in 1:length(Results$Compound)){
    ID = substr(Results_ratio$Compound[r],1,4)
    if(T %in% grepl(ID,IS$ID)){
      ISTD_ID = ID
    } else {
      ISTD = IS_Assignment[grepl(ID,IS_Assignment$Compound),]$ISTD
      ISTD_ID = substr(ISTD,1,4)
    }
    Values = as.numeric(Results[r,16:length(Results)])
    k = match(Results[grepl(ISTD_ID,Results$Compound),]$Compound,Results$Compound)
    Values_IS = as.numeric(Results[k,16:length(Results)])
    Ratio = Values / Values_IS
    Results_ratio[r,16:length(Results)] = Ratio
    if(Results_ratio$`Blank values`[r]>0){
      mean_Values_IS = mean(Values_IS,na.rm=T)
      if(!is.na(mean_Values_IS)|!is.nan(mean_Values_IS)){
        Results_ratio$`Blank values`[r] = Results_ratio$`Blank values`[r]/mean_Values_IS
      } else {
        Results_ratio$`Blank values`[r] = -1
      }
    }
  }
  
  #Test if there is an outlier in the references and remove it
  currentFunction = "Reference_Ratios"
  if(T%in%grepl(".mzML",files)){
    Reference_Ratios = Results_ratio %>% dplyr::select(contains(c("Compound","ISTD",str_remove(References,".mzML"))))
  } else {
    Reference_Ratios = Results_ratio %>% dplyr::select(contains(c("Compound","ISTD",str_remove(References,".raw"))))
  }
  
  if(T %in% grepl("CAL",References)){
    
    for(i in 1:length(Reference_Ratios$Compound)){
      Levels = gsub(".*CAL","CAL",colnames(Reference_Ratios)[3:length(Reference_Ratios)])
      names = Levels
      Levels = gsub("p",".",Levels)
      Levels_table = as.data.frame(table(Levels))
      Levels_table$Levels = names
      if(max(Levels_table$Freq)>1){
        for(j in 1:length(Levels_table$Levels)){
          relevant = Results_ratio %>% dplyr::select(contains(as.character(Levels_table$Levels[j])))
          quantiles = quantile(relevant[i,],na.rm = T)
          IQR = quantiles[4]-quantiles[2]
          low=quantiles[2]-1.5*IQR
          high=quantiles[4]+1.5*IQR
          Outlier_Reference = as.logical(relevant[i,]>high | relevant[i,]<low)
          if(is.na(Outlier_Reference)){
            Outlier_Reference = F
          }
          if(T%in%Outlier_Reference){
            Names = colnames(relevant)[Outlier_Reference]
            for(n in 1:length(Names)){
              k = match(Names[n],colnames(Reference_Ratios))
              Reference_Ratios[i,k] = NA
            }
          }
        }
      }
    }
  } else {
    PDF_path = paste0(results_path,"/",Sys.Date(),"_Reference_stability.pdf")
    pdf(file=PDF_path)
    for(i in 1:length(Reference_Ratios$Compound)){
      mean = mean(as.numeric(Reference_Ratios[i,3:length(Reference_Ratios)]),na.rm = T)
      sd = sd(as.numeric(Reference_Ratios[i,3:length(Reference_Ratios)]),na.rm = T)
      SEM = sd/sqrt(length(Reference_Ratios[i,3:length(Reference_Ratios)]))
      CV = sd/mean*100
      if(is.na(CV)|CV<=20){
        threshold = Reference_deviation
      } else {
        threshold = CV/100
      }
      low=mean-mean*threshold
      high=mean+mean*threshold
      #Use IQR if references are too variable (less affected by outliers)
      if(is.na(CV)|CV>(Reference_deviation*100)){
        quantiles = quantile(Reference_Ratios[i,3:length(Reference_Ratios)],na.rm = T)
        IQR = quantiles[4]-quantiles[2]
        low=quantiles[2]-0.75*IQR
        high=quantiles[4]+0.75*IQR
      }
      Outlier_Reference = as.logical(Reference_Ratios[i,3:length(Reference_Ratios)]>high | Reference_Ratios[i,3:length(Reference_Ratios)]<low)
      Outlier_Reference = Outlier_Reference[!is.na(Outlier_Reference)]
      if(length(Outlier_Reference)==0){
        Results_ratio$`Valid references`[i] = 0
        next
      }
      Results_ratio$`Valid references`[i] = length(Outlier_Reference[Outlier_Reference!=T])
      Reference_Ratios[i,c(rep(F,2),Outlier_Reference)] = NA
      Outlier_Reference_table = as.data.frame(table(Outlier_Reference))
      v = match(F,Outlier_Reference_table$Var1)
      if(T%in%Outlier_Reference){
        Results_ratio$Comment[i] = paste0(paste0(colnames(Reference_Ratios)[c(F,F,Outlier_Reference)],"_removed_as_outlier"),Results_ratio$Comment[i],collapse="|")
      }
      #plot
      if(Results_ratio$Compound[i]==Results_ratio$ISTD[i]){
        next
      }
      levels = as.numeric(as.numeric(Reference_Ratios[i,3:length(Reference_Ratios)]))
      if(mean(is.na(levels))==1){
        next
      }
      
      if(!is.na(CV)&CV<=20){
        deviation = 0.2*mean
      } else {
        deviation = 0.75*IQR
      }
      
      # Use CV as measure of variability
      plot_title <- paste(Reference_Ratios$Compound[i], "-", "CV:", round(CV, 2), "%")
      
      label = "Response Ratio"
      
      # Calculate the y-axis limits
      y_range <- range(levels)
      y_diff <- diff(y_range)
      deviation_limit <- ifelse(deviation > y_diff, 2 * deviation, 2 * y_diff)
      y_min <- min(levels) - deviation_limit
      y_max <- max(levels) + deviation_limit
      
      # Generate the plot
      p <- ggplot(data = NULL, aes(x = factor(seq_along(colnames(Reference_Ratios)[3:(length(Reference_Ratios))])), y = levels)) +
        geom_point(shape = 16, size = 4, color = "black") +
        geom_hline(yintercept = mean(levels,na.rm=T), linetype = "dashed", color = "red") +
        geom_hline(yintercept = mean(levels,na.rm = T) + deviation, linetype = "dashed", color = "orange") +
        geom_hline(yintercept = mean(levels,na.rm = T) - deviation, linetype = "dashed", color = "orange") +
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = low, ymax = high),
                  fill = "lightgreen", alpha = 0.2) +
        labs(x = NULL, y = label, title = plot_title) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold"),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
        scale_x_discrete(labels = colnames(Reference_Ratios)[3:(length(Reference_Ratios))]) +
        coord_cartesian(ylim = c(y_min, y_max))
      
      print(p)
    }
    dev.off()
  }
  
  for(i in 1:length(Results_ratio$Compound)){
    if(Results_ratio$Compound[i]==Results_ratio$ISTD[i]){
      Results_ratio$ISTD[i] = NA
    }
  }
  
  #Calculate Recovery (QA) or concentration (CAL)
  currentFunction = "calculate_recovery_or_concentration"
  wo_IS = Reference_Ratios[Reference_Ratios$Compound!=Reference_Ratios$ISTD,]
  Results_ratio = Results_ratio[!is.na(Results_ratio$ISTD),]
  
  if(T %in% grepl("CAL",References)){
    Results_ratio$R2_calibration = 1
    Results_ratio$`Valid calibration` = TRUE
    Levels = gsub(".*CAL","CAL",colnames(Reference_Ratios)[3:length(Reference_Ratios)])
    Levels = gsub("p",".",Levels)
    Conc = as.numeric(substr(Levels,5,8))
    models = list()
    length(models) = length(wo_IS$Compound)
    PDF_path = paste0(results_path,"/",Sys.Date(),"_Calibration_Curves.pdf")
    pdf(file=PDF_path,width=14,height=7)
    for(quan in 1:length(models)){
      count = 1
      modellas = list()
      r2s = 1
      accuracies = 1
      concentrations = 1
      influentials_list = list()
      Cutoffs = 1
      influentials = NA
      Cutoff = NA
      frame = data.frame("y"=Conc,"x"=as.numeric(wo_IS[quan,3:length(wo_IS)]))
      frame$x[frame$x==Inf|frame$x==-Inf] = 0
      frame = frame[frame$x!=0,]
      frame = frame[!is.na(frame$x),]
      if(length(frame$y[!is.na(frame$y)]) < 2){
        models[[quan]] = NA
        Results_ratio$R2_calibration[quan] = 0
        Results_ratio$`Valid calibration`[quan] = FALSE
        next
      }
      if(length(frame$y[frame$y!=0])<minimum_nr_of_datapoints_calibration){
        Results_ratio$Comment[quan] = paste0("Calibration points < 4",Results_ratio$Comment[quan],collapse = ";")
        Results_ratio$`Valid calibration`[quan] = FALSE
      }
      model = tryCatch(lm(formula=y~0+x,data = frame),error=function(e){model = NA})
      if(is.na(model)){
        models[[quan]] = NA
        Results_ratio$`Valid calibration`[quan] = FALSE
        Results_ratio$R2_calibration[quan] = 0
        next
      } 
      r2 = round(summary(model)$adj.r.squared,digits = 2)
      predi = predict(model,newdata = frame)
      accuracy = abs(median(((frame$y-predi)/frame$y)*100))
      redo_calib = F
      if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
        redo_calib = T
        r2_initial = r2
        accuracy_initial = accuracy
        while(redo_calib == T){
          if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
            break
          }
          cooksD = cooks.distance(model)
          influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          if(length(influential)==0){
            influential = NaN
          }
          if(!is.nan(influential)){
            frame_new = frame
            frame_new$y[match(max(influential),frame_new$y)] = NaN
            frame_new$y[frame_new$y%in%influential] = NaN
            frame_new = frame_new[!is.nan(frame_new$y),]
            model_new = lm(formula=y~0+x,data = frame_new)
            r2 = round(summary(model_new)$adj.r.squared,digits = 2)
            predi = predict(model_new,newdata = frame_new)
            accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
            if(r2<r2_initial&accuracy_initial<accuracy){
              redo_calib = F
              r2 = r2_initial
              accuracy = accuracy_initial
            } else {
              models[[quan]] = model_new
              if(r2 >= minimum_r2&accuracy<=30){
                redo_calib = F
              }
              frame = frame_new
              model = model_new
              r2_initial = r2
              accuracy_initial = accuracy
              influentials = append(influentials,influential)
            }
          } else {
            redo_calib = F
          }
        }
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      #using weights of x 
      ############################################################################
      model = tryCatch(lm(formula=y~0+x,data = frame,weights = x),error=function(e){model = NA})
      if(is.na(model)){
        models[[quan]] = NA
        Results_ratio$`Valid calibration`[quan] = FALSE
        Results_ratio$R2_calibration[quan] = 0
        next
      } 
      r2 = round(summary(model)$adj.r.squared,digits = 2)
      predi = predict(model,newdata = frame)
      accuracy = abs(median(((frame$y-predi)/frame$y)*100))
      redo_calib = F
      if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
        redo_calib = T
        r2_initial = r2
        accuracy_initial = accuracy
        while(redo_calib == T){
          if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
            break
          }
          cooksD = cooks.distance(model)
          influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          if(length(influential)==0){
            influential = NaN
          }
          if(!is.nan(influential)){
            frame_new = frame
            frame_new$y[match(max(influential),frame_new$y)] = NaN
            frame_new$y[frame_new$y%in%influential] = NaN
            frame_new = frame_new[!is.nan(frame_new$y),]
            model_new = lm(formula=y~0+x,data = frame_new,weights = x)
            r2 = round(summary(model_new)$adj.r.squared,digits = 2)
            predi = predict(model_new,newdata = frame_new)
            accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
            if(r2<r2_initial&accuracy_initial<accuracy){
              redo_calib = F
              r2 = r2_initial
              accuracy = accuracy_initial
            } else {
              models[[quan]] = model_new
              if(r2 >= minimum_r2&accuracy<=30){
                redo_calib = F
              }
              frame = frame_new
              model = model_new
              r2_initial = r2
              accuracy_initial = accuracy
              influentials = append(influentials,influential)
            }
          } else {
            redo_calib = F
          }
        }
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      
      #using weights of 1/x
      ############################################################################
      model = tryCatch(lm(formula=y~0+x,data = frame,weights = 1/x),error=function(e){model = NA})
      if(is.na(model)){
        models[[quan]] = NA
        Results_ratio$`Valid calibration`[quan] = FALSE
        Results_ratio$R2_calibration[quan] = 0
        next
      } 
      r2 = round(summary(model)$adj.r.squared,digits = 2)
      predi = predict(model,newdata = frame)
      accuracy = abs(median(((frame$y-predi)/frame$y)*100))
      redo_calib = F
      if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
        redo_calib = T
        r2_initial = r2
        accuracy_initial = accuracy
        while(redo_calib == T){
          if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
            break
          }
          cooksD = cooks.distance(model)
          influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          if(length(influential)==0){
            influential = NaN
          }
          if(!is.nan(influential)){
            frame_new = frame
            frame_new$y[match(max(influential),frame_new$y)] = NaN
            frame_new$y[frame_new$y%in%influential] = NaN
            frame_new = frame_new[!is.nan(frame_new$y),]
            model_new = lm(formula=y~0+x,data = frame_new,weights = 1/x)
            r2 = round(summary(model_new)$adj.r.squared,digits = 2)
            predi = predict(model_new,newdata = frame_new)
            accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
            if(r2<r2_initial&accuracy_initial<accuracy){
              redo_calib = F
              r2 = r2_initial
              accuracy = accuracy_initial
            } else {
              models[[quan]] = model_new
              if(r2 >= minimum_r2&accuracy<=30){
                redo_calib = F
              }
              frame = frame_new
              model = model_new
              r2_initial = r2
              accuracy_initial = accuracy
              influentials = append(influentials,influential)
            }
          } else {
            redo_calib = F
          }
        }
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      
      #without forcing through 0
      ############################################################################
      Cutoff = NA
      influentials = NA
      frame = data.frame("y"=Conc,"x"=as.numeric(wo_IS[quan,3:length(wo_IS)]))
      frame$x[frame$x==Inf|frame$x==-Inf] = 0
      frame = frame[frame$x!=0,]
      frame = frame[!is.na(frame$x),]
      model = tryCatch(lm(formula=y~x,data = frame,weights = x),error=function(e){model = NA})
      Cutoff = ((-model$coefficients[1])/model$coefficients[2])
      frame = frame[frame$x>=Cutoff,]
      model = tryCatch(lm(formula=y~x,data = frame),error=function(e){model = NA})
      
      if(!is.na(model)&length(frame$y)>=minimum_nr_of_datapoints_calibration){
        r2 = round(summary(model)$adj.r.squared,digits = 2)
        predi = predict(model,newdata = frame)
        accuracy = abs(median(((frame$y-predi)/frame$y)*100))
        redo_calib = F
        if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
          redo_calib = T
          r2_initial = r2
          accuracy_initial = accuracy
          while(redo_calib == T){
            if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
              break
            }
            cooksD = cooks.distance(model)
            influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
            if(length(influential)==0){
              influential = NaN
            }
            if(!is.nan(influential)){
              frame_new = frame
              frame_new$y[match(max(influential),frame_new$y)] = NaN
              frame_new = frame_new[!is.nan(frame_new$y),]
              model_new = lm(formula=y~x,data = frame_new)
              r2 = round(summary(model_new)$adj.r.squared,digits = 2)
              predi = predict(model_new,newdata = frame_new)
              accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
              if(r2<r2_initial|accuracy>accuracy_initial){
                redo_calib = F
                r2 = r2_initial
                accuracy = accuracy_initial
                influentials[influential%in%influentials] = NA
              } else {
                model = model_new
                if(r2 >= minimum_r2&accuracy<=30){
                  redo_calib = F
                }
                frame = frame_new
                model = model_new
                r2_initial = r2
                accuracy_initial = accuracy
                influentials = append(influentials,influential)
              }
            } else {
              redo_calib = F
            }
          }}
      } else {
        model = NA
        r2 = NA
        accuracy = NA
        Cutoff = NA
        influentials = NA
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      
      #without forcing through 0, weights = x
      ############################################################################
      Cutoff = NA
      influentials = NA
      frame = data.frame("y"=Conc,"x"=as.numeric(wo_IS[quan,3:length(wo_IS)]))
      frame$x[frame$x==Inf|frame$x==-Inf] = 0
      frame = frame[frame$x!=0,]
      frame = frame[!is.na(frame$x),]
      model = tryCatch(lm(formula=y~x,data = frame,weights=x),error=function(e){model = NA})
      Cutoff = ((-model$coefficients[1])/model$coefficients[2])
      frame = frame[frame$x>=Cutoff,]
      model = tryCatch(lm(formula=y~x,data = frame,weights = x),error=function(e){model = NA})
      
      if(!is.na(model)&length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
        r2 = round(summary(model)$adj.r.squared,digits = 2)
        predi = predict(model,newdata = frame)
        accuracy = abs(median(((frame$y-predi)/frame$y)*100))
        redo_calib = F
        if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
          redo_calib = T
          r2_initial = r2
          accuracy_initial = accuracy
          while(redo_calib == T){
            if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
              break
            }
            cooksD = cooks.distance(model)
            influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
            if(length(influential)==0){
              influential = NaN
            }
            if(!is.nan(influential)){
              frame_new = frame
              frame_new$y[match(max(influential),frame_new$y)] = NaN
              frame_new = frame_new[!is.nan(frame_new$y),]
              model_new = lm(formula=y~x,data = frame_new,weights = x)
              r2 = round(summary(model_new)$adj.r.squared,digits = 2)
              predi = predict(model_new,newdata = frame_new)
              accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
              if(r2<r2_initial&accuracy>accuracy_initial){
                redo_calib = F
                r2 = r2_initial
                accuracy = accuracy_initial
                influentials[influential%in%influentials] = NA
              } else {
                model = model_new
                if(r2 >= minimum_r2&accuracy<=30){
                  redo_calib = F
                }
                frame = frame_new
                r2_initial = r2
                accuracy_initial = accuracy
                influentials = append(influentials,influential)
              }
            } else {
              redo_calib = F
            }
          }}
      } else {
        model = NA
        r2 = NA
        accuracy = NA
        Cutoff = NA
        influentials = NA
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      
      #without forcing through 0, weights = 1/x
      ############################################################################
      Cutoff = NA
      influentials = NA
      frame = data.frame("y"=Conc,"x"=as.numeric(wo_IS[quan,3:length(wo_IS)]))
      frame$x[frame$x==Inf|frame$x==-Inf] = 0
      frame = frame[frame$x!=0,]
      frame = frame[!is.na(frame$x),]
      model = tryCatch(lm(formula=y~x,data = frame,weights=x),error=function(e){model = NA})
      Cutoff = ((-model$coefficients[1])/model$coefficients[2])
      frame = frame[frame$x>=Cutoff,]
      model = tryCatch(lm(formula=y~x,data = frame,weights = 1/x),error=function(e){model = NA})
      
      if(!is.na(model)&length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
        r2 = round(summary(model)$adj.r.squared,digits = 2)
        predi = predict(model,newdata = frame)
        accuracy = abs(median(((frame$y-predi)/frame$y)*100))
        redo_calib = F
        if((r2 < minimum_r2|accuracy>maximum_error_of_quantitation) & length(frame$y[frame$y!=0])>minimum_nr_of_datapoints_calibration){
          redo_calib = T
          r2_initial = r2
          accuracy_initial = accuracy
          while(redo_calib == T){
            if(length(frame$y[frame$y!=0])==minimum_nr_of_datapoints_calibration){
              break
            }
            cooksD = cooks.distance(model)
            influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
            if(length(influential)==0){
              influential = NaN
            }
            if(!is.nan(influential)){
              frame_new = frame
              frame_new$y[match(max(influential),frame_new$y)] = NaN
              frame_new = frame_new[!is.nan(frame_new$y),]
              model_new = lm(formula=y~x,data = frame_new,weights = 1/x)
              r2 = round(summary(model_new)$adj.r.squared,digits = 2)
              predi = predict(model_new,newdata = frame_new)
              accuracy = abs(median(((frame_new$y-predi)/frame_new$y)*100))
              if(r2<r2_initial&accuracy>accuracy_initial){
                redo_calib = F
                r2 = r2_initial
                accuracy = accuracy_initial
                influentials[influential%in%influentials] = NA
              } else {
                model = model_new
                if(r2 >= minimum_r2&accuracy<=30){
                  redo_calib = F
                }
                frame = frame_new
                r2_initial = r2
                accuracy_initial = accuracy
                influentials = append(influentials,influential)
              }
            } else {
              redo_calib = F
            }
          }}
      } else {
        model = NA
        r2 = NA
        accuracy = NA
        Cutoff = NA
        influentials = NA
      }
      modellas[[count]] = model
      r2s[count] = r2
      if(!is.na(model)){
        concentrations[count] = length(summary(model)$residuals)
      } else {
        concentrations[count] = 0
      }
      accuracies[count] = accuracy
      Cutoffs[count] = Cutoff
      influentials_list[[count]] = influentials
      count = count+1
      
      indices = which(accuracies<maximum_error_of_quantitation&r2s>minimum_r2)
      indices2 = indices[which(concentrations[indices]==max(concentrations[indices]))]
      if(length(indices2)>1){
        aki = match(min(accuracies[indices2],na.rm=T),accuracies)
      } else {
        aki = match(accuracies[indices2],accuracies)
      }
      if(length(aki)==0){
        aki = match(min(accuracies,na.rm=T),accuracies)
      }
      accuracy = accuracies[aki]
      model = modellas[[aki]]
      r2 = r2s[aki]
      Cutoff = Cutoffs[aki]
      if(length(model$coefficients)>1){
        Cutoff_y = min(model$model$y[model$model$y>=model$coefficients[1]&model$model$y>0])
      } else {
        Cutoff_y = min(model$model$y[model$model$y>0])
      }
      Cutoff_y_max = max(model$model$y)
      influentials = influentials_list[[aki]]
      models[[quan]] = model
      
      if(!is.nan(r2)){
        if(is.na(model)){
          models[[quan]] = NA
          Results_ratio$`Valid calibration`[quan] = FALSE
          Results_ratio$R2_calibration[quan] = 0
          next
        }
        
        if(r2<minimum_r2|accuracy>maximum_error_of_quantitation){
          Results_ratio$Comment[quan] = paste0("r2 < ",minimum_r2," and/or accuracy > ",maximum_error_of_quantitation," %",Results_ratio$Comment[quan],collapse = ";")
          Results_ratio$R2_calibration[quan] = r2
          Results_ratio$`Valid calibration`[quan] = FALSE
          Results_ratio$Min_Calibration[quan] = Cutoff_y
          Results_ratio$Max_Calibration[quan] = max(model$model$y)
        } else {
          Results_ratio$R2_calibration[quan] = r2
          Results_ratio$Min_Calibration[quan] = Cutoff_y
          Results_ratio$Max_Calibration[quan] = max(model$model$y)
          #Results_ratio$`Valid calibration`[quan] = TRUE
        }
      }
      mylabel = bquote(italic(R)^2 == .(format(r2, digits = 2)))
      frame = data.frame("y"=Conc,"x"=as.numeric(wo_IS[quan,3:length(wo_IS)]))
      frame$x[frame$x==Inf|frame$x==-Inf] = 0
      frame = frame[frame$x!=0,]
      frame = frame[!is.na(frame$x),]
      frame$col = rep("black",length(frame$y))
      if(length(influentials[!is.na(influentials)])>0){
        frame$col[frame$y%in%influentials[!is.na(influentials)]] = "red"
      }
      if(!is.na(Cutoff)){
        frame$col[frame$x<=Cutoff] = "red"
        frame$col[frame$y<Cutoff_y] = "red"
      }
      
      data <- data.frame(x = frame$x, y = frame$y, col = frame$col)
      
      title <- str_wrap(paste0("Compound: ", wo_IS$Compound[quan], " - Standard: ", wo_IS$ISTD[quan]), width = 14*72)
      
      p = ggplot(data, aes(x = x, y = y, colour = col)) +
        geom_point(size=3) +                       
        labs(x = "Response Ratio", y = paste0("Conc. [",quantified_unit,"]"),
             title = title) +
        geom_abline(intercept = ifelse(length(coef(model)) > 1, coef(model)[1], 0),
                    slope = coef(model)[length(coef(model))],
                    color = "blue") +
        xlim(c(-0.005, max(frame$x, na.rm = TRUE) * 1.2)) +
        ylim(c(0, max(frame$y, na.rm = TRUE) * 1.2)) +
        geom_vline(xintercept = Cutoff, col = "red") +
        geom_hline(yintercept = Cutoff_y, col = "red") +
        geom_hline(yintercept = Cutoff_y_max, col = "red") +
        annotate("text", x = 0, y = max(frame$y, na.rm = TRUE) * 1.2,
                 label = mylabel, size = 5, hjust = 0, vjust = 1) +
        theme_minimal() +   # Use a minimalistic theme
        scale_color_manual(values = c("black"="black","red"="red")) +
        theme(axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 18),
              legend.position = "none",
              plot.margin = margin(t = 20, unit = "pt"))
      print(p)
      
    }
    dev.off()
    Results_final = Results_ratio
    Results_final$IS = NULL
    
    for(quant in 1:length(Results_final$Compound)){
      frame = data.frame("x"=as.numeric(Results_final[quant,16:length(Results_final)]))
      if(is.na(models[[quant]])){
        Results_final[quant,16:length(Results_final)] = NA
        next
      }
      quantified = predict(models[[quant]],newdata = frame)
      Results_final[quant,16:length(Results_final)] = quantified
      if(Results_final$`Blank values`[quant]>0){
        if(!is.na(model)){
          frame = data.frame("x"=as.numeric(Results_final$`Blank values`[quant]))
          Results_final$`Blank values`[quant]  = as.numeric(predict(models[[quant]],newdata = frame))
        } else {
          Results_final$`Blank values`[quant]  = -1
        }
      }
      
    }
  } else {
    Results_final = Results_ratio
    for(rec in 1:length(Results_final$Compound)){
      mean = as.numeric(rowMeans(wo_IS[rec,3:length(wo_IS)],na.rm = T))
      Results_final[rec,16:length(Results_final)] = (as.numeric(Results_final[rec,16:length(Results_final)])/mean)*100
      
      if(Results_final$`Blank values`[rec]>0&!is.na(mean)&!is.nan(mean)){
        Results_final$`Blank values`[rec] = as.numeric(Results_final$`Blank values`[rec])/mean*100
      } else if (Results_final$`Blank values`>0) {
        Results_final$`Blank values`[rec] = -1
      }
    }
    
  }
  
  #Add the problematic peaks -> Get "CHECK" in excluded
  for(a in 16:ncol(Results_final)){
    Peaklist_a = Peaklists[[a-15]]
    for(b in 1:nrow(Results_final)){
      c = match(Results_final$Compound[b],Peaklist_a$Compound)
      d = match(Results_final$Compound[b],Results$Compound)
      if(grepl("<a7>",Peaklist_a$Comment[c])){
        area_a = stringr::str_extract(Peaklist_a$Comment[c],"<a7>*<a7>")
        area_a = as.numeric(str_match(Peaklist_a$Comment[c], "<a7>\\s*(.*?)\\s*<a7>")[,2])
        IS_a = Results$ISTD[c]
        k = match(IS_a,Results$Compound)
        Area_IS_a = as.numeric(Results[k,a])
        ratio_a = data.frame("x"=area_a/Area_IS_a)
        if(!grepl("CAL",References_refined)){
          mean = as.numeric(rowMeans(wo_IS[b,3:length(wo_IS)],na.rm = T))
          Results_final[b,a] = (as.numeric(ratio_a)/mean)*100
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
        } else {
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
          if(!is.na(models[[b]])){
            value_a = as.numeric(predict(models[[b]],ratio_a))
            Results_final[b,a] = value_a
          }
        }
      }
    }
  }
  
  #Add the qualitative peaks -> Get "QUAL" in excluded
  currentFunction = "add_qualitative_or_problematic_peaks_to_peaks_all"
  
  for(a in 16:ncol(Results_final)){
    Peaklist_a = Peaklists[[a-15]]
    for(b in 1:nrow(Results_final)){
      c = match(Results_final$Compound[b],Peaklist_a$Compound)
      d = match(Results_final$Compound[b],Results$Compound)
      if(grepl("?",Peaklist_a$Comment[c],fixed = T)){
        area_a = stringr::str_extract(Peaklist_a$Comment[c],"\\?*\\?")
        area_a = as.numeric(str_match(Peaklist_a$Comment[c], "\\?\\s*(.*?)\\s*\\?")[,2])
        IS_a = Results$ISTD[c]
        k = match(IS_a,Results$Compound)
        Area_IS_a = as.numeric(Results[k,a])
        ratio_a = data.frame("x"=area_a/Area_IS_a)
        if(!grepl("CAL",References_refined)){
          mean = as.numeric(rowMeans(wo_IS[b,3:length(wo_IS)],na.rm = T))
          Results_final[b,a] = (as.numeric(ratio_a)/mean)*100
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
        } else {
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
          if(!is.na(models[[b]])){
            value_a = as.numeric(predict(models[[b]],ratio_a))
            Results_final[b,a] = value_a
          }
        }
      }
    }
  }
  
  #Add potential negative peaks -> Dont get character in excluded
  for(a in 16:ncol(Results_final)){
    Peaklist_a = Peaklists[[a-15]]
    for(b in 1:nrow(Results_final)){
      c = match(Results_final$Compound[b],Peaklist_a$Compound)
      d = match(Results_final$Compound[b],Results$Compound)
      if(grepl("&",Peaklist_a$Comment[c])){
        area_a = stringr::str_extract(Peaklist_a$Comment[c],"&*&")
        area_a = as.numeric(str_match(Peaklist_a$Comment[c], "&\\s*(.*?)\\s*&")[,2])
        IS_a = Results$ISTD[c]
        k = match(IS_a,Results$Compound)
        Area_IS_a = as.numeric(Results[k,a])
        ratio_a = data.frame("x"=area_a/Area_IS_a)
        if(!grepl("CAL",References_refined)){
          mean = as.numeric(rowMeans(wo_IS[b,3:length(wo_IS)],na.rm = T))
          Results_final[b,a] = (as.numeric(ratio_a)/mean)*100
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
        } else {
          Results[d,a] = as.numeric(area_a)
          Results_ratio[b,a] = as.numeric(ratio_a)
          if(!is.na(models[[b]])){
            value_a = as.numeric(predict(models[[b]],ratio_a))
            Results_final[b,a] = value_a
          }
        }
      }
    }
  }
  
  Output = openxlsx::createWorkbook()
  if(use.area==T){
    openxlsx::addWorksheet(Output,"Areas")
    openxlsx::writeData(Output,"Areas",Results)
  } else {
    openxlsx::addWorksheet(Output,"Heights")
    openxlsx::writeData(Output,"Heights",Results)
  }
  
  openxlsx::addWorksheet(Output,"Ratios")
  openxlsx::writeData(Output,"Ratios",Results_ratio)
  openxlsx::addWorksheet(Output,"Final")
  openxlsx::writeData(Output,"Final",Results_final)
  
  Results_final_all = Results_final
  
  currentFunction = "create_Results_final_excluded"
  
  if(T %in% grepl("CAL",References)){
    for(n in 1:length(Results_final$Compound)){
      if(!is.na(models[[n]])&!is.na(Results_final$R2_calibration[n])&Results_final$R2_calibration[n]!=0){
        maxi = max(models[[n]]$model$x)
        y_intercept = models[[n]]$coefficients[1]
        y_min = Results_final$Min_Calibration[n]
        slope = models[[n]]$coefficients[2]
        mini = models[[n]]$model$x[match(y_min,models[[n]]$model$y)]
        Ratios = Results_ratio[n,16:length(Results_final)]
        Concis = Results_final[n,16:length(Results_final)]
        NA_cases = is.na(Ratios)
        Ratios[is.na(Ratios)] = 0
        Concis[Ratios > (maxi*max_calibration)] = "above calibration"
        Concis[Ratios <= mini] = "below calibration"
        Concis[NA_cases] = NA
        Row = cbind(Results_final[n,(1:15)],Concis)
        Results_final[n,] = Row
      }
    }
  } else {
    names = colnames(Results_final)
    Reference_files = names[grepl(Reference_ID,names)]
    for(n in 1:nrow(Results_final)){
      excluded = as.logical(sapply(Reference_files,grepl,x=Results_final$Comment[n]))
      for(m in 1:length(Reference_files)){
        if(excluded[m]==T){
          rut = match(Reference_files[m],names)
          Results_final[n,rut] = "Removed"
        }
      }
    }
  }
  
  for(a in 16:ncol(Results_final)){
    Peaklist_a = Peaklists[[a-15]]
    for(b in 1:nrow(Results_final)){
      c = match(Results_final$Compound[b],Peaklist_a$Compound)
      if(grepl("<a7>",Peaklist_a$Comment[c],fixed=T)&!grepl(Solvent_Blank_ID,colnames(Results_final)[a])){
        Results_final[b,a] = "CHECK"
      }
    }
  }
  
  if(!(T%in%grepl("QA",References_refined))){
    for(a in 16:ncol(Results_final)){
      Peaklist_a = Peaklists[[a-15]]
      for(b in 1:nrow(Results_final)){
        c = match(Results_final$Compound[b],Peaklist_a$Compound)
        if(grepl("?",Peaklist_a$Comment[c],fixed = T)){
          Results_final[b,a] = "QUAL"
        }
      }
    }
  }
  
  
  for(b in 1:length(Results_final$Compound)){
    if(Results_final$`Blank values`[b]!=0){
      ha = match(Results_final$Compound[b],Results_qualitative$Compound)
      Filter = as.logical(Results_qualitative[ha,2:length(Results_qualitative)]=="Background")
      Filter[is.na(Filter)] = F
      Results_final[b,c(rep(F,15),Filter)] = "masked by background"
    }
  }
  
  openxlsx::addWorksheet(Output,"Final_excluded")
  openxlsx::writeData(Output,"Final_excluded",Results_final)
  
  openxlsx::addWorksheet(Output,"Qualitative")
  openxlsx::writeData(Output,"Qualitative",Results_qualitative)
  TRUE_style = openxlsx::createStyle(fontColour = "#000000", bgFill = "lightgreen")
  FALSE_style = openxlsx::createStyle(fontColour = "#000000", bgFill = "lightcoral")
  Background_style = openxlsx::createStyle(fontColour = "#000000", bgFill = "orange")
  openxlsx::conditionalFormatting(Output,"Qualitative", cols = 1:ncol(Results_qualitative),
                                  rows = 1:nrow(Results_qualitative)+1,rule = "FALSE",style = FALSE_style,
                                  type = "contains")
  openxlsx::conditionalFormatting(Output,"Qualitative", cols = 1:ncol(Results_qualitative),
                                  rows = 1:nrow(Results_qualitative)+1,rule = "TRUE",style = TRUE_style,
                                  type = "contains")
  openxlsx::conditionalFormatting(Output,"Qualitative", cols = 1:ncol(Results_qualitative),
                                  rows = 1:nrow(Results_qualitative)+1,rule = "Background",style = Background_style,
                                  type = "contains")
  openxlsx::addWorksheet(Output,"Peaklist")
  openxlsx::writeData(Output,"Peaklist",Peaks)
  openxlsx::addWorksheet(Output,"Retentiontimes")
  openxlsx::writeData(Output,"Retentiontimes",Results_RT)
  openxlsx::saveWorkbook(Output,paste0(results_path,"/",Sys.Date(),"_",sample,"_","Results_ATS.xlsx"),overwrite = T)
  
  if(length(Intensity_dependent_shift)>0){
    IS_shift_wb = openxlsx::createWorkbook()
    for(nun in 1:length(Intensity_dependent_shift)){
      name = names(Intensity_dependent_shift)[nun]
      name = stringr::str_replace(name,":",".")
      if(nchar(name)>31){
        name = stringr::str_sub(name,1,31)
      }
      openxlsx::addWorksheet(IS_shift_wb,name)
      openxlsx::writeData(IS_shift_wb,name,Intensity_dependent_shift[[nun]])
    }
    openxlsx::saveWorkbook(IS_shift_wb,paste0(results_path,"/",Sys.Date(),"_",sample,"_","Intensity_dependent_shifts.xlsx"),overwrite = T)
  }
  
  combinations = as.data.frame(tidyr::crossing(Results_final$Compound,colnames(Results_final)[16:ncol(Results_final)]))
  colnames(combinations) = c("Compound","File")
  combinations$Potential_value = NA
  for(zui in 1:nrow(combinations)){
    a = match(combinations$Compound[zui],Results_final_all$Compound)
    b = match(combinations$File[zui],colnames(Results_final_all))
    combinations$Potential_value[zui] = Results_final_all[a,b]
  }
  combinations$Current_value = NA
  for(zui in 1:nrow(combinations)){
    a = match(combinations$Compound[zui],Results_final$Compound)
    b = match(combinations$File[zui],colnames(Results_final))
    combinations$Current_value[zui] = Results_final[a,b]
  }
  combinations$Action = NA
  combinations$Start_RT = NA
  combinations$End_RT = NA
  
  for_manual = openxlsx::createWorkbook()
  openxlsx::addWorksheet(for_manual,"What to do")
  openxlsx::addWorksheet(for_manual,"worksheet")
  Actions = c("NA/Empty",
              "Accept",
              "Remove",
              "Redo")
  Reaction = c("nothing happens",
               "current_value will be replaced by potential_value",
               "current_value will be set to NA",
               "you have to set a Start_ and End_RT and it will get the Area and calculate respective outcomes")
  instructions = data.frame(Actions,Reaction)
  openxlsx::writeData(for_manual,"What to do",instructions)
  openxlsx::writeData(for_manual,"worksheet",combinations)
  openxlsx::saveWorkbook(for_manual,paste0(results_path,"/",Sys.Date(),"_",sample,"_for_manual_evaluation.xlsx"),overwrite = T)
  
  dir.create(paste0(results_path,"/results_plots"))
  
  currentFunction = "writing PDFs"
  print("writing PDFs")
  pb <- txtProgressBar(max=nrow(Results), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  pdf_frame = matrix(seq(1,nrow(Results),1),ncol=cores)
  
  rest = nrow(Results)-nrow(pdf_frame)*cores
  
  if(rest!=0){
    for(i in (nrow(pdf_frame)+rest):nrow(pdf_frame)){
      pdf_frame[i,ncol(pdf_frame)] = NA
    }
  }
  
  if(cores==1){
    for(p in 1:nrow(Results)){
      if(T%in%(Peaks$Compound[p]%in%worksheet_new$Compound)){
        qpdf::pdf_combine(input=mixedsort(list.files(paste0(results_path,"/plots"),full.names = T,pattern = gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[p])))),
                          output=paste0(results_path,"/results_plots/",gsub("[()]","",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[p]))),".pdf"))
      }
    }
  } else {
    pb <- txtProgressBar(max=cores, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    cl <- makeCluster(length(References))
    doSNOW::registerDoSNOW(cl)
    environment(errors) = environment()
    tryCatch({
      foreach(i=1:cores,
              .packages = c("qpdf","stringr","gtools"),
              .options.snow=opts,
              .export=c("results_path","Peaks")) %dopar% {
                for(p in pdf_frame[,i]){
                  if(is.na(p)){
                    next
                  }
                  qpdf::pdf_combine(input=mixedsort(list.files(paste0(results_path,"/plots"),full.names = T,pattern = gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[p])))),
                                    output=paste0(results_path,"/results_plots/",gsub("[()]","",gsub("\\[|\\]"," ",gsub("[()]"," ",Peaks$Compound[p]))),".pdf"))
                }
              }
      close(pb)
      stopCluster(cl)
    },error=errors)
  }
  save.image(file=paste0(results_path,"/",Sys.Date(),"_R_workspace.RData"))
  
  Output_function = list()
  Output_function$Results_ratio = Results_ratio
  Output_function$Results_final_all =Results_final_all
  Output_function$Results_final = Results_final
  return(Output_function)
}

