###############################################################################
#                                                                        Feb 19
#          Discrete Choice Stan Model - RBT Prey Selection
#
#  Notes:
#  
#  To Do:
#  * reduce duplication in the if else cases 
#
###############################################################################




model.set.up = function(model.name = NULL){
  #-----------------------------------------------------------------------------#
  # read in the drift data
  
  if(is.null(model.name)){
    message("Need a model name...")
  } else {
    tmp.A = read.table(paste0(getwd(), "/Data/Drift_A.txt"), sep = "\t", header = TRUE)
    colnames(tmp.A)[3:80] = as.character(1:78) 
    
    # only by taxa, not size
    tmp.A.w = read.table(paste0(getwd(), "/Data/Drift_A_w.txt"), sep = "\t", header = TRUE)
    
    # these should match
    if(identical(tmp.A[,1:2], tmp.A.w[,1:2]) == FALSE){
      message("Stop - Problem with the data")
    }
    #-----------------------------------------------------------------------------#
    # read in the diet data
    
    y.tmp2 = read.table(paste0(getwd(), "/Data/Diet.txt"), sep = "\t", header = TRUE)
    colnames(y.tmp2)[6:83] = as.character(1:78) 
    
    # only by taxa, not size
    w.tmp2 = read.table(paste0(getwd(), "/Data/Diet_w.txt"), sep = "\t", header = TRUE)
    
    # these should match
    if(identical(y.tmp2[,1:3], w.tmp2[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    } 
    #-----------------------------------------------------------------------------#
    # format/match/check data
    
    drift.ts = paste(tmp.A[,1], tmp.A[,2])
    diet.ts = sort(unique(paste(y.tmp2[,1], y.tmp2[,2])))
    
    # subset only the drift that matches the diet
    A = as.matrix(tmp.A[match(diet.ts, drift.ts),3:ncol(tmp.A)])
    w.a = as.matrix(tmp.A.w[match(diet.ts, drift.ts),3:ncol(tmp.A.w)])
    Nst = as.numeric(nrow(A))  # number of site & trips
    
    # diet
    y.tmp = y.tmp2[order(paste(y.tmp2[,1], y.tmp2[,2])),]
    w.tmp = w.tmp2[order(paste(w.tmp2[,1], w.tmp2[,2])),]
    
    y.tmp <<- y.tmp # writes this to the global
    
    if(identical(y.tmp[,1:3], w.tmp[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    }
    y.in = as.matrix(y.tmp[,6:ncol(y.tmp)])
    w.in = as.matrix(w.tmp[,4:ncol(w.tmp)])
    
    #-----------------------------------------------------------------------------#
    # format some variables for the model
    Nind = nrow(y.in)   # Number of individuals
    Nsp = 7             # number of prey taxa 
    upper = c(7, 15, 10, 13, 8, 12, 20) - 1
    Nspsz = sum(upper)  # number of taxa and sizes
    spsz = c(1:Nspsz)
    idx_ts_ind = as.numeric(as.factor(paste(y.tmp$trip, y.tmp$site)))
    
    spp = sort(as.numeric(rep(c(1:Nsp), upper)))
    
    # new varible for the fixed paramaters in the drift portion of the model 
    spp2 = sort(as.numeric(rep(c(1:Nsp), upper-1)))
    
    # better way to do this? 
    idx = c(rep(1, upper[1]),
            rep(upper[1] + 1, upper[2]),
            rep(upper[1] + upper[2] + 1, upper[3]),
            rep(upper[1] + upper[2] + upper[3] + 1, upper[4]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + 1, upper[5]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + 1, upper[6]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6] + 1, upper[7]))
    
    idx2 = c(rep(upper[1], upper[1]),
             rep(upper[1] + upper[2], upper[2]),
             rep(upper[1] + upper[2] + upper[3], upper[3]),
             rep(upper[1] + upper[2] + upper[3] + upper[4], upper[4]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5], upper[5]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6], upper[6]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6] + upper[7], upper[7]))
    
    out = NULL
    for(i in 1:Nsp){
      out[[i]] = seq(2, upper[i] + 1)  
    }
    
    sz = do.call(c,out)
    
    len.idx = sz
    len.idx <<- len.idx
    
    u.idx = unique(idx)
    u.idx2 = unique(idx2)
    
    # need these 
    Ns = 5
    Nt = length(unique(y.tmp[,1]))
    trip = as.numeric(y.tmp[,1])
    site = as.numeric(y.tmp[,2])
    
    #-----------------------------------------------------------------------------#
    # fish size covariate (one sample has no size, set to the mean)
    t.sz = y.tmp[,4] 
    
    t.sz[which(t.sz == 2287)] = NA
    
    fish_sz = t.sz - mean(t.sz, na.rm = T) 
    
    fish_sz[is.na(fish_sz)] = 0
    
    #-----------------------------------------------------------------------------#
    # need emp_a (drift prop for each taxa across all the data)
    tmp.a = data.frame(spp = spp,
                       sums = colSums(A))
    
    tmp.a.2 = group_by(tmp.a, spp) %>%
      summarize(my.sum = sum(sums))
    
    a.tmp.w = colSums(w.a)
    
    tmp.a.all = tmp.a.2$my.sum + a.tmp.w
    
    emp_a_2 = tmp.a.all / sum(tmp.a.all)
    
    # emp.dat2 = data.frame(taxa = name.key$name,
    #                       emp = emp_a_2)
    #-----------------------------------------------------------------------------#
    if(model.name == "Length"){
      # calc the average size of prey in the drift (probably a better way :/)
      temp = data.frame(sz = sz,
                        sums = colSums(A))
      
      temp.2 = group_by(temp, sz) %>%
        summarize(my.sum = sum(sums))
      
      tmp.vec = list()
      
      for(i in 1:dim(temp.2)[1]){
        tmp.vec[[i]] = rep(as.character(temp.2[i,1]), times = as.numeric(temp.2[i,2]))
      }
      
      all.len = do.call('c', tmp.vec)
      
      avg.log.measure = mean(log(as.numeric(all.len)))
    }
    
    #-----------------------------------------------------------------------------#
    if(model.name == "Width"){
      # width calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput' on 'out.w'
      
      fit.dat.w = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                                 V2 = c(-0.712264417633861, -1.585825466801, -1.33393917278193,
                                        -1.80644981811059, -1.34684364955096, -0.295195222715158, 
                                        -1.49260182871226),
                                 V3 = c(0.589397273790521, 0.800389362070142, 0.614460012015383,
                                        0.563246021013028, 0.890778931760663, 0.424033122560842, 
                                        0.137416988547813)),
                            row.names = c(NA, -7L), class = "data.frame")
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat.w[which(fit.dat.w$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  width = exp(tmp.ar))
      }
      
      dat.width = do.call('rbind', area.df)
      
      avg.log.measure = mean(log(dat.width$width))
      measure = dat.width$width
    }
    #-----------------------------------------------------------------------------#
    if(model.name == "Area"){
      # area calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput'
      
      fit.dat = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                               V2 = c(0.224460916479209, -1.78662907514783, -1.34743919150295,
                                      -2.15070453700794, -1.45525542402561, -0.74863545618833,
                                      -1.85570798696119),
                               V3 = c(0.94958327948254, 1.79612553310148, 1.54243803018606,
                                      1.67798706868551, 1.89261925347815, 1.65868421232674,
                                      1.30725930407769)),
                          row.names = c(NA, -7L), class = "data.frame")
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat[which(fit.dat$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  area = exp(tmp.ar))
      }
      
      # area = do.call('c', tmp.ar2)
      area = do.call('rbind', area.df)
      avg.log.measure = mean(log(area$area))
      measure = area$area
      
    }    
    #-----------------------------------------------------------------------------#
    if(model.name == "Mass"){
      # estimate mass (see foodbase package species list for regression parameters)
      # species list to estimate biomass
      sp.tmp = read.table(paste0(getwd(), "/Data/SpeciesList_2019_02_07.txt"), sep = "\t", header = TRUE)
      
      sp.tmp.2 = sp.tmp[which(sp.tmp$SpeciesID %in% c("CHIL", "CHIA", "SIML", "SIMA",
                                                      "GAMM", "NZMS", "OLIG")),]
      sp.tmp.3 = sp.tmp.2[,which(names(sp.tmp.2) %in% c("SpeciesID", "RegressionA", "RegressionB"))]
      
      sp.key = data.frame(num = c(2,4,5,7,19,21,1),
                          name = c("NZMS", "GAMM", "SIMA", "SIML", "CHIA", "CHIL", "OLIG"))
      sp.key$A = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionA
      sp.key$B = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionB
      
      
      out.list = list() 
      
      for(i in 1:Nsp){
        tmp.sz = sz[which(spp == i)]
        
        out.list[[i]] = sp.key[i,"A"] * tmp.sz^sp.key[i,"B"]
        
      }
      
      # mass
      measure = do.call(c,out.list)
      avg.log.measure = mean(log(measure))
      
      # check / look at regressions
      # test = as.data.frame(cbind(sz, as.factor(spp), mass))
      # windows()
      # p = ggplot(test, aes(x = sz, y = mass, group = spp)) +
      #     geom_point(aes(color = spp)) + 
      #     geom_line(aes(color = spp)) 
      # p
      
      
      #-----------------------------------------------------------------------------#  
    }
    
    # dirichlet params 
    alpha = rep(.2, Nsp)
    
    # dummy matrix
    X <- as.matrix(model.matrix(~ as.factor(spp) - 1))   
    
    idx_first = unique(idx)  # position of the first bin for each taxa
    tmpper = seq(1:Nspsz)
    not_first = tmpper[which(!tmpper %in% idx_first)]
    #-----------------------------------------------------------------------------#
    # add this in above 
    if(model.name == "Length"){
      measure = sz
    }
    
    data.in = list(Nspsz = Nspsz, Nst = Nst, Nsp = Nsp, Nind = Nind,
                   sp = spp, sp_idx2 = spp2, spsz = spsz,
                   a = A, w_a = w.a, y = y.in, w = w.in,
                   idx = idx, idx2 = idx2, idx_first = idx_first, not_first = not_first,
                   idx_ts_ind = idx_ts_ind,
                   alpha = alpha, a_Nsz = upper,   
                   X = X, 
                   emp_a = emp_a_2,
                   sz = measure,  u_idx = u.idx, u_idx2 = u.idx2,
                   avg_log_len = avg.log.measure, 
                   fish_sz = fish_sz)  
    
    return(data.in)  
  }
}

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
model.set.up.no.worms = function(model.name = NULL){
  #-----------------------------------------------------------------------------#
  # read in the drift data
  
  if(is.null(model.name)){
    message("Need a model name...")
  } else {
    tmp.A = read.table(paste0(getwd(), "/Data/Drift_A.txt"), sep = "\t", header = TRUE)
    colnames(tmp.A)[3:80] = as.character(1:78) 
    
    # cut out worms
    tmp.A = tmp.A[,-c((ncol(tmp.A) - 18):ncol(tmp.A))]
    
    # only by taxa, not size
    tmp.A.w = read.table(paste0(getwd(), "/Data/Drift_A_w.txt"), sep = "\t", header = TRUE)
    
    # cut out worms
    tmp.A.w = tmp.A.w[,-9]
    
    # these should match
    if(identical(tmp.A[,1:2], tmp.A.w[,1:2]) == FALSE){
      message("Stop - Problem with the data")
    }
    #-----------------------------------------------------------------------------#
    # read in the diet data
    
    y.tmp2 = read.table(paste0(getwd(), "/Data/Diet.txt"), sep = "\t", header = TRUE)
    colnames(y.tmp2)[6:83] = as.character(1:78) 
    
    # cut out the worms
    y.tmp2 = y.tmp2[,-c(60:78)]
    
    # only by taxa, not size
    w.tmp2 = read.table(paste0(getwd(), "/Data/Diet_w.txt"), sep = "\t", header = TRUE)
    
    # cut out the worms
    w.tmp2 = w.tmp2[,-10]
    
    # these should match
    if(identical(y.tmp2[,1:3], w.tmp2[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    } 
    
    #-----------------------------------------------------------------------------#
    # format/match/check data
    
    drift.ts = paste(tmp.A[,1], tmp.A[,2])
    diet.ts = sort(unique(paste(y.tmp2[,1], y.tmp2[,2])))
    
    # subset only the drift that matches the diet
    A = as.matrix(tmp.A[match(diet.ts, drift.ts),3:ncol(tmp.A)])
    w.a = as.matrix(tmp.A.w[match(diet.ts, drift.ts),3:ncol(tmp.A.w)])
    Nst = as.numeric(nrow(A))  # number of site & trips
    
    # diet
    y.tmp = y.tmp2[order(paste(y.tmp2[,1], y.tmp2[,2])),]
    w.tmp = w.tmp2[order(paste(w.tmp2[,1], w.tmp2[,2])),]
    
    y.tmp <<- y.tmp  # writes this to the global
    
    if(identical(y.tmp[,1:3], w.tmp[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    }
    y.in = as.matrix(y.tmp[,6:ncol(y.tmp)])
    w.in = as.matrix(w.tmp[,4:ncol(w.tmp)])
    
    #-----------------------------------------------------------------------------#
    # format some variables for the model
    Nind = nrow(y.in)   # Number of individuals
    Nsp = 6             # number of prey taxa 
    upper = c(7, 15, 10, 13, 8, 12) - 1
    Nspsz = sum(upper)  # number of taxa and sizes
    spsz = c(1:Nspsz)
    idx_ts_ind = as.numeric(as.factor(paste(y.tmp$trip, y.tmp$site)))
    
    spp = sort(as.numeric(rep(c(1:Nsp), upper)))
    
    # new varible for the fixed paramaters in the drift portion of the model 
    spp2 = sort(as.numeric(rep(c(1:Nsp), upper-1)))
    
    # better way to do this? 
    idx = c(rep(1, upper[1]),
            rep(upper[1] + 1, upper[2]),
            rep(upper[1] + upper[2] + 1, upper[3]),
            rep(upper[1] + upper[2] + upper[3] + 1, upper[4]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + 1, upper[5]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + 1, upper[6]))
    
    idx2 = c(rep(upper[1], upper[1]),
             rep(upper[1] + upper[2], upper[2]),
             rep(upper[1] + upper[2] + upper[3], upper[3]),
             rep(upper[1] + upper[2] + upper[3] + upper[4], upper[4]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5], upper[5]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6], upper[6]))
    
    out = NULL
    for(i in 1:Nsp){
      out[[i]] = seq(2, upper[i] + 1)  
    }
    
    sz = do.call(c,out)
    
    len.idx = sz
    len.idx <<- len.idx
    
    u.idx = unique(idx)
    u.idx2 = unique(idx2)
    
    # need these (not for the current model)
    # Ns = 5
    # Nt = length(unique(y.tmp[,1]))
    # trip = as.numeric(y.tmp[,1])
    # site = as.numeric(y.tmp[,2])
    
    #-----------------------------------------------------------------------------#
    # fish size covariate (one sample has no size, set to the mean)
    t.sz = y.tmp[,4] 
    
    t.sz[which(t.sz == 2287)] = NA
    
    fish_sz = t.sz - mean(t.sz, na.rm = T) 
    
    fish_sz[is.na(fish_sz)] = 0
    
    #-----------------------------------------------------------------------------#
    # need emp_a (drift prop for each taxa across all the data)
    tmp.a = data.frame(spp = spp,
                       sums = colSums(A))
    
    tmp.a.2 = group_by(tmp.a, spp) %>%
      summarize(my.sum = sum(sums))
    
    a.tmp.w = colSums(w.a)
    
    tmp.a.all = tmp.a.2$my.sum + a.tmp.w
    
    emp_a_2 = tmp.a.all / sum(tmp.a.all)
    
    # emp.dat2 = data.frame(taxa = name.key$name,
    #                       emp = emp_a_2)
    #-----------------------------------------------------------------------------#
    if(model.name == "Length"){
      # calc the average size of prey in the drift (probably a better way :/)
      temp = data.frame(sz = sz,
                        sums = colSums(A))
      
      temp.2 = group_by(temp, sz) %>%
        summarize(my.sum = sum(sums))
      
      tmp.vec = list()
      
      for(i in 1:dim(temp.2)[1]){
        tmp.vec[[i]] = rep(as.character(temp.2[i,1]), times = as.numeric(temp.2[i,2]))
      }
      
      all.len = do.call('c', tmp.vec)
      
      avg.log.measure = mean(log(as.numeric(all.len)))
    }
    
    #-----------------------------------------------------------------------------#
    if(model.name == "Width"){
      # width calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput' on 'out.w'
      
      fit.dat.w = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                                 V2 = c(-0.712264417633861, -1.585825466801, -1.33393917278193,
                                        -1.80644981811059, -1.34684364955096, -0.295195222715158, 
                                        -1.49260182871226),
                                 V3 = c(0.589397273790521, 0.800389362070142, 0.614460012015383,
                                        0.563246021013028, 0.890778931760663, 0.424033122560842, 
                                        0.137416988547813)),
                            row.names = c(NA, -7L), class = "data.frame")
      # cut out worms
      fit.dat.w = fit.dat.w[-7,]
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      name.key = name.key[-7,]
      
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat.w[which(fit.dat.w$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  width = exp(tmp.ar))
      }
      
      dat.width = do.call('rbind', area.df)
      
      avg.log.measure = mean(log(dat.width$width))
      measure = dat.width$width
    }
    #-----------------------------------------------------------------------------#
    if(model.name == "Area"){
      # area calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput'
      
      fit.dat = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                               V2 = c(0.224460916479209, -1.78662907514783, -1.34743919150295,
                                      -2.15070453700794, -1.45525542402561, -0.74863545618833,
                                      -1.85570798696119),
                               V3 = c(0.94958327948254, 1.79612553310148, 1.54243803018606,
                                      1.67798706868551, 1.89261925347815, 1.65868421232674,
                                      1.30725930407769)),
                          row.names = c(NA, -7L), class = "data.frame")
      
      fit.dat = fit.dat[-7,]
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      name.key = name.key[-7,]
      
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat[which(fit.dat$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  area = exp(tmp.ar))
      }
      
      # area = do.call('c', tmp.ar2)
      area = do.call('rbind', area.df)
      avg.log.measure = mean(log(area$area))
      measure = area$area
      
    }    
    #-----------------------------------------------------------------------------#
    if(model.name == "Mass"){
      # estimate mass (see foodbase package species list for regression parameters)
      # species list to estimate biomass
      sp.tmp = read.table(paste0(getwd(), "/Data/SpeciesList_2019_02_07.txt"), sep = "\t", header = TRUE)
      
      sp.tmp.2 = sp.tmp[which(sp.tmp$SpeciesID %in% c("CHIL", "CHIA", "SIML", "SIMA",
                                                      "GAMM", "NZMS")),]
      sp.tmp.3 = sp.tmp.2[,which(names(sp.tmp.2) %in% c("SpeciesID", "RegressionA", "RegressionB"))]
      
      sp.key = data.frame(num = c(2,4,5,7,19,21),
                          name = c("NZMS", "GAMM", "SIMA", "SIML", "CHIA", "CHIL"))
      sp.key$A = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionA
      sp.key$B = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionB
      
      out.list = list() 
      
      for(i in 1:Nsp){
        tmp.sz = sz[which(spp == i)]
        
        out.list[[i]] = sp.key[i,"A"] * tmp.sz^sp.key[i,"B"]
        
      }
      
      # mass
      measure = do.call(c,out.list)
      avg.log.measure = mean(log(measure))
      
      # check / look at regressions
      # test = as.data.frame(cbind(sz, as.factor(spp), mass))
      # windows()
      # p = ggplot(test, aes(x = sz, y = mass, group = spp)) +
      #     geom_point(aes(color = spp)) + 
      #     geom_line(aes(color = spp)) 
      # p
      
      
      #-----------------------------------------------------------------------------#  
    }
    
    # dirichlet params 
    alpha = rep(.2, Nsp)
    
    # dummy matrix
    X <- as.matrix(model.matrix(~ as.factor(spp) - 1))   
    
    idx_first = unique(idx)  # position of the first bin for each taxa
    tmpper = seq(1:Nspsz)
    not_first = tmpper[which(!tmpper %in% idx_first)]
    #-----------------------------------------------------------------------------#
    # add this in above 
    if(model.name == "Length"){
      measure = sz
    }
    
    data.in = list(Nspsz = Nspsz, Nst = Nst, Nsp = Nsp, Nind = Nind,
                   sp = spp, sp_idx2 = spp2, spsz = spsz,
                   a = A, w_a = w.a, y = y.in, w = w.in,
                   idx = idx, idx2 = idx2, idx_first = idx_first, not_first = not_first,
                   idx_ts_ind = idx_ts_ind,
                   alpha = alpha, a_Nsz = upper,   
                   X = X, 
                   emp_a = emp_a_2,
                   sz = measure,  u_idx = u.idx, u_idx2 = u.idx2,
                   avg_log_len = avg.log.measure, 
                   fish_sz = fish_sz)  
    
    return(data.in)  
  }
}
#-----------------------------------------------------------------------------#
