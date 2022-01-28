#' Calculate Haplotypes Differences between population
#'
#' @export calculate_haplotypic_differences
#' @param filename .arp file intended to be used for the analysis
#'
#' @examples create_haplo_pop_dataset("example.arp")
#' @examples create_haplo_pop_dataset("~/dir/example.arp")
#'
#' @return A text summary file and diagnostic plot gets saved in the current directory for the model.
#' @return A density plot is also saved in the environment to show the Haplotype differences per region.
#' @import reshape2 ggplot2 ggfortify stringr
#' @description This package builds a model to test if populations and regions can explain the differences in haplotypes sequences between individuals. You need to provide a .arp file with subject from different regions and their haplotypes sequence

###R Only
#Haplo to data frame
#Read the file by line

calculate_haplotypic_differences <- function(filename){

    ##create vector with all the line in the file
    vec1 <- readLines(filename)
    #vector with our tag selection
    parts <- c("\tHaplList={", "[[Samples]]", "[[Structure]] ")

    #Boolean vector, return TRUE at the line that in v1 that contain 1 element of our vector
    vec2 <-vec1 %in% parts

    #create dataframe with all the line and all the result
    df1 <- data.frame(vec1,vec2)
    #transform line in a list of character
    df1$vec1 <- as.character(df1$vec1)
    #cumulative sum of our vec2, on each TRUE = incremente by 1
    df1$vec3 <- cumsum(df1$vec2)


    #################################################################################
    #create the haplo dataset with haplotype data
    haplo_dataset <- df1[df1$vec3 ==1,]
    parts <- c("\tHaplList={","}")
    haplo_dataset$vec2 <-haplo_dataset$vec1 %in% parts
    haplo_dataset$vec3 <-cumsum(haplo_dataset$vec2)
    haplo_dataset <- haplo_dataset[haplo_dataset$vec3 == 1,]
    haplo_dataset <- as.data.frame(as.character(haplo_dataset[-1,1]))
    colnames(haplo_dataset) <- "Haplotype"
    haplo_dataset <- cbind(colsplit(haplo_dataset$Haplotype, pattern = "  ",names = c("ID", " Haplotype")))
    haplo_dataset$ID <- formatC(haplo_dataset$ID,width=4,format="d", flag=0)
    haplo_dataset$ID <- as.factor(haplo_dataset$ID)

    # to save the output dataframe in R, use the <<- to assign it to a dataframe, so we save the haplotypes dataset
    haplotypes <- haplo_dataset

    ####################
    #population dataset
    grp_name_abrev <- df1[df1$vec3 == 3,]


    pop_dataset <- df1[df1$vec3 == 2,]



    #Group_name vectore that will be used to have the abreviation
    group_name <- pop_dataset$vec1[startsWith(pop_dataset$vec1, "#\t")]

    # remove the #/t from the beginning of each group name
    group_name <- lapply(group_name, function (x){
      gsub("#\t", "", x)
    })

    # create the abbreviation of each group name
    AAA_gname <- c()
    for (gname in 1:length(group_name)) {
      A <- substr(strsplit(unlist(group_name)[gname], split = " ")[[1]][1],start = 1,stop = 1)
      AA <- substr(strsplit(unlist(group_name)[gname], split = " ")[[1]][2],start = 1,stop = 2)
      AAA_gname <- append(x= AAA_gname, values = paste(A,AA, sep = ""), after = gname)
    }

    ##Sample_name vector
    sample_name <- pop_dataset$vec1[startsWith(pop_dataset$vec1, "\tSampleName")]
    #clean the vector
    sample_name <- lapply(sample_name,function (x) {
      gsub("\tSampleName= \"", "", x)
    })

    sample_name <- lapply(sample_name, function(x){
      gsub("\"","",x)
    })

    #Sample data vector
    sample_data <- pop_dataset$vec1
    sample_data[startsWith(sample_data, "#\t")] <- TRUE
    sample_data <- sample_data[grepl("[[:digit:]]", sample_data) | grepl(TRUE, sample_data)]


    #remove samplesize frome the vector
    sample_data <- sample_data[!grepl( "\tSampleSize",sample_data)]


    #Create a vector where all the element that contain # REF will have TRUE
    vec <- c(grepl("# REF", sample_data))

    #Transform the sample_data vector to a dataframe
    sample_data <- as.data.frame(sample_data)

    #add the TRUE / FALSE vector to our dataframe, that will serv as a base
    sample_data$Ref_tag <- vec

    #Create a vector that will match our sample_data vector lenght and will contain the pop_name
    #with the exact occurence between the TRUE.
    pop_name_match_data <- c(NA)
    pop_index <- 0
    for(i in sample_data$Ref_tag){
      if(i == TRUE){

        pop_index <- pop_index + 1

      }
      pop_name_match_data <- append(pop_name_match_data, sample_name[pop_index], after = length(pop_name_match_data))
    }
    pop_name_match_data <- as.character(pop_name_match_data)



    ##Same principle to add the abreviation:
    pop_abreviation <- c()
    abrev_index <- 0
    for(i in sample_data$sample_data){
      if (i == TRUE){
        abrev_index <- abrev_index+1
      }
      pop_abreviation <- append(pop_abreviation, AAA_gname[abrev_index], after = length(pop_abreviation))
    }

    #Add it to the dataframe and clean the dataframe from the TRUE value (# REF rows)
    sample_data$pop_name <- pop_name_match_data
    sample_data$Group_abrev <- pop_abreviation
    sample_data<- sample_data[!(sample_data$Ref_tag == TRUE),]
    sample_data <- sample_data[!(sample_data$sample_data == TRUE),]
    sample_data$pop_name <- factor(sample_data$pop_name)
    sample_data$Group_abrev <- as.factor(sample_data$Group_abrev)
    sample_data <- sample_data[c(1,3,4)]

    sample_data_final <-  cbind(colsplit(sample_data$sample_data,  "\     ",names = c("ID", "Sample_Size")))
    sample_data_final$Pop_name <- sample_data$pop_name
    sample_data_final$region <- sample_data$Group_abrev


    sample_data_final$ID <- formatC(sample_data_final$ID,width=4,format="d", flag=0)
    sample_data_final<- sample_data_final[c(2,3,4,1)]
    # View(sample_data_final)

    ##Dup the sample_size to get pop_id (2 haplo 0002 -> 1:0002, 2:0002)
    #new data for safety and simplicity
    final_dataset_test <- data.frame(Pop_id = numeric(),
                                     Pop_name = character(),
                                     region = factor(),
                                     ID = factor())

    index <- 1

    for (i in sample_data_final$Sample_Size){
      y <- 1
      final_dataset_test <- rbind(final_dataset_test, sample_data_final[index,])
      if(i > 1){
        while (y < i){
          final_dataset_test <- rbind(final_dataset_test, sample_data_final[index,])

          y <- y + 1
        }
      }
      index <- index +1

    }

    #rename samplesize to pop_id
    colnames(final_dataset_test)<- c("Pop_id","Pop_name","region", "ID")

    #ID_pop
    sample_ID <- c()

    for (number in 1:sum(as.numeric(sample_data_final$Sample_Size))) {
      sample_ID <- append(sample_ID, number)
    }

    final_dataset_test$Pop_id <- sample_ID

    colnames(final_dataset_test) <- c("individual_id", "pop_name", "region","haplo_id")

    # to save the output dataframe in R, use the <<- to assign it to a dataframe, so we save the whole dataset
    final_dataset <- final_dataset_test


  countdifferences <- function(haplotypes){
    difference <- data.frame(matrix(0,nrow(haplotypes),nrow(haplotypes)))
    for (i in 1:nrow(haplotypes)){
      for (j in 2:nrow(haplotypes)){
        seq1 = unlist(strsplit(haplotypes[i,2],""))
        seq2 = unlist(strsplit(haplotypes[j,2],""))
        if (j>i){ #to avoid to compare the same haplotype and to be sure to have only once each pair of haplotype
          diff = sum(seq1 != seq2)
          difference[i,j] = diff
          difference[j,i] = diff
        }
      }
    }
    names(difference) = c(haplotypes[,1]) #sso that the names of the columns are the same as the ref of the haplotypes
    rownames(difference) = c(haplotypes[,1])
    return(difference)
  }

  HaploDiff <- countdifferences(haplotypes)


  pop <- final_dataset
  pop <- as.data.frame(paste(unlist(pop$individual_id), pop$pop_name, pop$region, pop$haplo_id, sep = "_"))

  combi <- expand.grid(pop[,1],pop[,1])
  ID1 <- as.data.frame(str_split_fixed(combi[,1], "_",n=4))[1]
  ID2 <- as.data.frame(str_split_fixed(combi[,2], "_", n=4))[1]
  combiunique <- combi[ID1 > ID2,] #Sp that we have only unique combination

  pop1 <- as.data.frame(str_split_fixed(combiunique$Var1, "_",n=4))
  pop2 <- as.data.frame(str_split_fixed(combiunique$Var2, "_", n=4))

  popcomp <- data.frame(matrix(0,nrow(pop1),3))
  popcomp[pop1[,2] != pop2[,2],1] <- ("Different_populations")
  popcomp[pop1[,2] == pop2[,2],1] <- pop1[pop1[,2] == pop2[,2],2]
  popcomp[pop1[,3] != pop2[,3],2] <- ("Different_regions")
  popcomp[pop1[,3] == pop2[,3],2] <- pop1[pop1[,3] == pop2[,3],3]
  colnames(popcomp) <- c("popnames", "regionnames", "haplotypes") #sso that the names of the columns are the same as the ref of the haplotypes

  popdiff <- function(i){
    HaploDiff[pop1[i,4],pop2[i,4]] #haplotypes differences taken from the matrix created using the function "countdifferences"
  }
  popcomp$haplotypes <- lapply(c(1:nrow(popcomp)), popdiff)

  # change the population and region names into factors, and unlist the haplotypes column
  popcomp$haplotypes <- unlist(popcomp$haplotypes)
  popcomp$popnames <- as.factor(popcomp$popnames)
  popcomp$regionnames <- as.factor(popcomp$regionnames)
  haplotypic_differences_df <<- popcomp

  # make model using lm(), print summary, then save it as a txt file
  linear_model <<- lm(data=popcomp, haplotypes ~ popnames + regionnames)
  sink("linear-model_summary.txt")
  print(summary(linear_model))
  sink()
  print(summary(linear_model))

  # make autoplot for model and save it in current directory
  png(filename = "model_diagnostic_plots.png", width = 1000, height = 719)
  autoplot <- autoplot(linear_model, which = 1:6)
  print(autoplot)
  dev.off()

  # make density plot for haplotypes by region names and save it in current directory
  ggplot(popcomp, aes(x = haplotypes, color = regionnames)) +
    geom_density(bw=0.9)+
    labs(title = "Density Plot for Haplotype Differences per Region",x="Number of Haplotype Differences",y="Density",color="Region Names") +
    theme_classic()
  ggsave("haplotypes_density_per_region.png")
}

