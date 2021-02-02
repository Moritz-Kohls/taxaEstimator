#' Class specific positive predictive values (PPV)
#'
#' Step 1 of the taxa classification and estimation procedure.
#' This function creates artificial reads based on NCBI viral reference genomes FASTA file,
#' computes input features and artificial neural network (ANN) model
#' and finally stores the trained model and classification results in R-intern variables.
#'
#' @seealso \code{\link{ANN_new_sample}}, \code{\link{estimation_taxa_distribution}}
#' @import Biostrings
#' @import keras
#' @import stringi
#' @param fasta.file_path FASTA file path (e.g. FASTA file viral.genomic.fna downloaded from NCBI: ftp://ftp.ncbi.nih.gov/refseq/release/viral).
#' @param taxonomy.file_path Taxonomy file path (e.g. taxonomy file taxonomy_viruses_available.csv delivered within this package).
#' @param temp.directory Results directory containing accuracy results, generalised confusion matrix results and accuracy as well as loss graphics of the simulation runs.
#' @param count.reads_training Number of sampled viruses and artificially generated reads per virus taxonomy, e.g. order (training data).
#' @param read_length Read length of all artificially generated reads.
#' @param simulation_runs Simulation runs
#' @examples
#' # Please specify your file paths and directories!
#' fasta.file_path = "~/ag_bioinf/genomes/viruses_na/refseq/viral.genomic.fna" # Download from NCBI!
#' taxonomy.file_path = "inst/extdata/taxonomy_viruses_available.csv" # Relative file path
#' temp.directory = "~/ag_bioinf/research/metagenomics/temp" # Results directory
#' count.reads_training = 100
#' read_length = 150
#' simulation_runs = 10
#' \dontrun{
#' class_specific_PPVs ( fasta.file_path, taxonomy.file_path, temp.directory, count.reads_training, read_length )
#' }
class_specific_PPVs = function ( fasta.file_path, taxonomy.file_path, temp.directory, count.reads_training, read_length, simulation_runs ) {

  # Taxonomy file of available viruses:
  taxonomy.file = read.csv(taxonomy.file_path, sep=";", stringsAsFactors=TRUE)
  str(taxonomy.file)
  taxonomy.file$NC_Id = as.character(taxonomy.file$NC_Id)
  taxonomy.file$Species = as.character(taxonomy.file$Species)
  str(taxonomy.file)
  table(taxonomy.file$Order, useNA = "always")

  # Taxonomy file of viruses with available order information:
  df.taxonomy_order = taxonomy.file[is.na(taxonomy.file$Order) == F,]

  # Fasta file:
  fasta.file = readDNAStringSet(fasta.file_path, format = "fasta")
  head(names(fasta.file))
  genome_lengths = width(fasta.file)
  head(genome_lengths)
  indices = regexpr(" ",names(fasta.file))
  fasta.NC_Ids = substring(names(fasta.file),1,indices-1)

  count.viruses_all = table(df.taxonomy_order$Order) # All viruses of the 9 different orders
  count.viruses_order = length(count.viruses_all) # 9 orders
  count.viruses_validation = count.viruses_test = round ( count.viruses_all * 0.15 ) # Number of viruses used for validation resp. test data
  count.viruses_training = count.viruses_all - count.viruses_validation - count.viruses_test # Number of viruses used for training data

  count.reads_validation = count.reads_test = round ( count.reads_training / 70 * 15 ) # Number of sampled viruses and artificially generated reads per virus taxonomy, e.g. order (validation resp. test data)
  count.reads_all = count.reads_training + count.reads_validation + count.reads_test # Number of sampled viruses and artificially generated reads per virus taxonomy, e.g. order (training, validation and test data)

  accuracy_results = rep(NA_real_,simulation_runs)
  confusion_matrix_results = vector("list",simulation_runs)
  acc_loss_graphics_results = vector("list",simulation_runs)

  for ( iteration in 1:(simulation_runs) ) {

    print(iteration)
    set.seed(iteration)
    # Split all viruses of all orders into training, validation and test data:
    Order_NC_Id.all = vector("list", length = count.viruses_order)
    for ( i in 1:count.viruses_order ) {
      Order_NC_Id.all [[i]] = unname(unlist(subset(df.taxonomy_order, Order == names(count.viruses_all)[i], select = NC_Id)))
    }
    names(Order_NC_Id.all) = names(count.viruses_all)
    Order_NC_Id.training = Order_NC_Id.validation = Order_NC_Id.test = vector("list", length = count.viruses_order)
    for ( i in 1:count.viruses_order ) {
      indices.training = sample(count.viruses_all[i], count.viruses_training[i])
      indices.temp = setdiff(1:count.viruses_all[i], indices.training)
      indices.validation = sample(indices.temp, count.viruses_validation[i])
      indices.test = setdiff(indices.temp, indices.validation)
      Order_NC_Id.training [[i]] = Order_NC_Id.all[[i]] [indices.training]
      Order_NC_Id.validation [[i]] = Order_NC_Id.all[[i]] [indices.validation]
      Order_NC_Id.test [[i]] = Order_NC_Id.all[[i]] [indices.test]
    }
    names(Order_NC_Id.training) = names(Order_NC_Id.validation) = names(Order_NC_Id.test) = names(count.viruses_all)

    # NC Ids of viruses per order. From each order, count.reads_training, count.reads_validation and count.reads_test NC Ids are sampled.
    order_NC_Id_sampled.training = matrix(nrow = count.viruses_order, ncol = count.reads_training)
    order_NC_Id_sampled.validation = matrix(nrow = count.viruses_order, ncol = count.reads_validation)
    order_NC_Id_sampled.test = matrix(nrow = count.viruses_order, ncol = count.reads_test)
    rownames(order_NC_Id_sampled.training) = rownames(order_NC_Id_sampled.validation) = rownames(order_NC_Id_sampled.test) = names(count.viruses_all)
    for ( i in 1:count.viruses_order ) {
      if ( count.viruses_training[i] >= count.reads_training ) {
        order_NC_Id_sampled.training [i,] = sample(Order_NC_Id.training [[i]], size = count.reads_training)
      }
      if ( count.viruses_validation[i] >= count.reads_validation ) {
        order_NC_Id_sampled.validation [i,] = sample(Order_NC_Id.validation [[i]], size = count.reads_validation)
      }
      if ( count.viruses_test[i] >= count.reads_test ) {
        order_NC_Id_sampled.test [i,] = sample(Order_NC_Id.test [[i]], size = count.reads_test)
      }
      if ( count.viruses_training[i] < count.reads_training ) {
        NC_Ids = rep(Order_NC_Id.training [[i]], floor(count.reads_training/count.viruses_training[i]))
        {
          if ( count.reads_training %% count.viruses_training[i] > 0 ) {
            order_NC_Id_sampled.training [i,] = c(NC_Ids, sample(Order_NC_Id.training [[i]], size = count.reads_training %% count.viruses_training[i]))
          }
          else if ( count.reads_training %% count.viruses_training[i] == 0 ) {
            order_NC_Id_sampled.training [i,] = NC_Ids
          }
        }
        order_NC_Id_sampled.training [i,] = sample(order_NC_Id_sampled.training [i,])
      }
      if ( count.viruses_validation[i] < count.reads_validation ) {
        NC_Ids = rep(Order_NC_Id.validation [[i]], floor(count.reads_validation/count.viruses_validation[i]))
        {
          if ( count.reads_validation %% count.viruses_validation[i] > 0 ) {
            order_NC_Id_sampled.validation [i,] = c(NC_Ids, sample(Order_NC_Id.validation [[i]], size = count.reads_validation %% count.viruses_validation[i]))
          }
          else if ( count.reads_validation %% count.viruses_validation[i] == 0 ) {
            order_NC_Id_sampled.validation [i,] = NC_Ids
          }
        }
        order_NC_Id_sampled.validation [i,] = sample(order_NC_Id_sampled.validation [i,])
      }
      if ( count.viruses_test[i] < count.reads_test ) {
        NC_Ids = rep(Order_NC_Id.test [[i]], floor(count.reads_test/count.viruses_test[i]))
        {
          if ( count.reads_test %% count.viruses_test[i] > 0 ) {
            order_NC_Id_sampled.test [i,] = c(NC_Ids, sample(Order_NC_Id.test [[i]], size = count.reads_test %% count.viruses_test[i]))
          }
          else if ( count.reads_test %% count.viruses_test[i] == 0 ) {
            order_NC_Id_sampled.test [i,] = NC_Ids
          }
        }
        order_NC_Id_sampled.test [i,] = sample(order_NC_Id_sampled.test [i,])
      }
    }

    # Indices (positions) of sampled NC Ids in fasta file:
    indices_NC_Id.training = matrix(sapply(order_NC_Id_sampled.training, FUN = function(x) grep(x,fasta.NC_Ids)), nrow = count.viruses_order, ncol = count.reads_training)
    indices_NC_Id.validation = matrix(sapply(order_NC_Id_sampled.validation, FUN = function(x) grep(x,fasta.NC_Ids)), nrow = count.viruses_order, ncol = count.reads_validation)
    indices_NC_Id.test = matrix(sapply(order_NC_Id_sampled.test, FUN = function(x) grep(x,fasta.NC_Ids)), nrow = count.viruses_order, ncol = count.reads_test)

    # Random start positions of reads:
    art_reads.start_pos.training = matrix(nrow = count.viruses_order, ncol = count.reads_training)
    art_reads.start_pos.validation = matrix(nrow = count.viruses_order, ncol = count.reads_validation)
    art_reads.start_pos.test = matrix(nrow = count.viruses_order, ncol = count.reads_test)
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_training ) {
        art_reads.start_pos.training [i,j] = sample.int(genome_lengths [ indices_NC_Id.training[i,j] ] - read_length + 1, size = 1)
      }
    }
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_validation ) {
        art_reads.start_pos.validation [i,j] = sample.int(genome_lengths [ indices_NC_Id.validation[i,j] ] - read_length + 1, size = 1)
      }
    }
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_test ) {
        art_reads.start_pos.test [i,j] = sample.int(genome_lengths [ indices_NC_Id.test[i,j] ] - read_length + 1, size = 1)
      }
    }

    # Generate artificial reads by subsequencing original reads:
    art_reads.training = vector("list", length = count.viruses_order * count.reads_training)
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_training ) {
        art_reads.training[[count.reads_training*(i-1)+j]] = subseq(fasta.file[[indices_NC_Id.training[i,j]]], start = art_reads.start_pos.training[i,j], width = read_length)
      }
    }
    art_reads.validation = vector("list", length = count.viruses_order * count.reads_validation)
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_validation ) {
        art_reads.validation[[count.reads_validation*(i-1)+j]] = subseq(fasta.file[[indices_NC_Id.validation[i,j]]], start = art_reads.start_pos.validation[i,j], width = read_length)
      }
    }
    art_reads.test = vector("list", length = count.viruses_order * count.reads_test)
    for ( i in 1:count.viruses_order ) {
      for ( j in 1:count.reads_test ) {
        art_reads.test[[count.reads_test*(i-1)+j]] = subseq(fasta.file[[indices_NC_Id.test[i,j]]], start = art_reads.start_pos.test[i,j], width = read_length)
      }
    }

    # One-mer, two-mer and three-mer distributions and inter-nucleotide distances (4 + 16 + 64 + 36 = 120 variables).
    two_mer_permutations = expand.grid(c("A","C","G","T"),c("A","C","G","T"))
    two_mer_permutations = as.matrix(two_mer_permutations)
    two_mer_permutations = apply(two_mer_permutations, MARGIN = 1, FUN = function (x) paste(x,collapse=""))
    two_mer_permutations = sort(two_mer_permutations)
    three_mer_permutations = expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T"))
    three_mer_permutations = as.matrix(three_mer_permutations)
    three_mer_permutations = apply(three_mer_permutations, MARGIN = 1, FUN = function (x) paste(x,collapse=""))
    three_mer_permutations = sort(three_mer_permutations)

    count_input_features = 4 + 16 + 64 + 36

    x_training = matrix(nrow = count.viruses_order * count.reads_training, ncol = count_input_features)
    for ( i in 1:(count.viruses_order*count.reads_training)) {
      x_training [i,1:4] = letterFrequency(art_reads.training[[i]], letters = c("A","C","G","T"), as.prob = F)
      for ( j in 5:20 ) {
        x_training [i,j] = stri_count_fixed(art_reads.training[[i]],pattern = two_mer_permutations[j-4])
      }
      for ( j in 21:84 ) {
        x_training [i,j] = stri_count_fixed(art_reads.training[[i]],pattern = three_mer_permutations[j-20])
      }
      dist_A = diff(stri_locate_all(art_reads.training[[i]], regex = "A") [[1]] [,1])
      dist_A = table(factor(dist_A, levels = 2:10))
      dist_C = diff(stri_locate_all(art_reads.training[[i]], regex = "C") [[1]] [,1])
      dist_C = table(factor(dist_C, levels = 2:10))
      dist_G = diff(stri_locate_all(art_reads.training[[i]], regex = "G") [[1]] [,1])
      dist_G = table(factor(dist_G, levels = 2:10))
      dist_T = diff(stri_locate_all(art_reads.training[[i]], regex = "T") [[1]] [,1])
      dist_T = table(factor(dist_T, levels = 2:10))
      x_training [i,85:93] = dist_A
      x_training [i,94:102] = dist_C
      x_training [i,103:111] = dist_G
      x_training [i,112:120] = dist_T
    }
    for ( my.index in 85:count_input_features ) {
      x_training[is.na(x_training[,my.index]),my.index] = 0
    }
    x_training = apply(x_training, MARGIN = 2, FUN = function(x) (x-min(x)) / diff(range(x)))
    colnames(x_training) = c("A","C","G","T",two_mer_permutations,three_mer_permutations,
                             paste0("d_A_",2:10),paste0("d_C_",2:10),paste0("d_G_",2:10),paste0("d_T_",2:10))

    x_validation = matrix(nrow = count.viruses_order * count.reads_validation, ncol = count_input_features)
    for ( i in 1:(count.viruses_order*count.reads_validation)) {
      x_validation [i,1:4] = letterFrequency(art_reads.validation[[i]], letters = c("A","C","G","T"), as.prob = F)
      for ( j in 5:20 ) {
        x_validation [i,j] = stri_count_fixed(art_reads.validation[[i]],pattern = two_mer_permutations[j-4])
      }
      for ( j in 21:84 ) {
        x_validation [i,j] = stri_count_fixed(art_reads.validation[[i]],pattern = three_mer_permutations[j-20])
      }
      dist_A = diff(stri_locate_all(art_reads.validation[[i]], regex = "A") [[1]] [,1])
      dist_A = table(factor(dist_A, levels = 2:10))
      dist_C = diff(stri_locate_all(art_reads.validation[[i]], regex = "C") [[1]] [,1])
      dist_C = table(factor(dist_C, levels = 2:10))
      dist_G = diff(stri_locate_all(art_reads.validation[[i]], regex = "G") [[1]] [,1])
      dist_G = table(factor(dist_G, levels = 2:10))
      dist_T = diff(stri_locate_all(art_reads.validation[[i]], regex = "T") [[1]] [,1])
      dist_T = table(factor(dist_T, levels = 2:10))
      x_validation [i,85:93] = dist_A
      x_validation [i,94:102] = dist_C
      x_validation [i,103:111] = dist_G
      x_validation [i,112:120] = dist_T
    }
    for ( my.index in 85:count_input_features ) {
      x_validation[is.na(x_validation[,my.index]),my.index] = 0
    }
    x_validation = apply(x_validation, MARGIN = 2, FUN = function(x) (x-min(x)) / diff(range(x)))
    colnames(x_validation) = c("A","C","G","T",two_mer_permutations,three_mer_permutations,
                               paste0("d_A_",2:10),paste0("d_C_",2:10),paste0("d_G_",2:10),paste0("d_T_",2:10))

    x_test = matrix(nrow = count.viruses_order * count.reads_test, ncol = count_input_features)
    for ( i in 1:(count.viruses_order*count.reads_test)) {
      x_test [i,1:4] = letterFrequency(art_reads.test[[i]], letters = c("A","C","G","T"), as.prob = F)
      for ( j in 5:20 ) {
        x_test [i,j] = stri_count_fixed(art_reads.test[[i]],pattern = two_mer_permutations[j-4])
      }
      for ( j in 21:84 ) {
        x_test [i,j] = stri_count_fixed(art_reads.test[[i]],pattern = three_mer_permutations[j-20])
      }
      dist_A = diff(stri_locate_all(art_reads.test[[i]], regex = "A") [[1]] [,1])
      dist_A = table(factor(dist_A, levels = 2:10))
      dist_C = diff(stri_locate_all(art_reads.test[[i]], regex = "C") [[1]] [,1])
      dist_C = table(factor(dist_C, levels = 2:10))
      dist_G = diff(stri_locate_all(art_reads.test[[i]], regex = "G") [[1]] [,1])
      dist_G = table(factor(dist_G, levels = 2:10))
      dist_T = diff(stri_locate_all(art_reads.test[[i]], regex = "T") [[1]] [,1])
      dist_T = table(factor(dist_T, levels = 2:10))
      x_test [i,85:93] = dist_A
      x_test [i,94:102] = dist_C
      x_test [i,103:111] = dist_G
      x_test [i,112:120] = dist_T
    }
    for ( my.index in 85:count_input_features ) {
      x_test[is.na(x_test[,my.index]),my.index] = 0
    }
    x_test = apply(x_test, MARGIN = 2, FUN = function(x) (x-min(x)) / diff(range(x)))
    colnames(x_test) = c("A","C","G","T",two_mer_permutations,three_mer_permutations,
                         paste0("d_A_",2:10),paste0("d_C_",2:10),paste0("d_G_",2:10),paste0("d_T_",2:10))

    # Categories to predict (Virus orders 1 to 9, here indices 0 to 8):
    y_training = rep(0:(count.viruses_order-1),each = count.reads_training)
    y_validation = rep(0:(count.viruses_order-1),each = count.reads_validation)
    y_test = rep(0:(count.viruses_order-1),each = count.reads_test)
    df.order_id = data.frame(Order = names(count.viruses_all), Id = 0:(count.viruses_order-1))

    # Build the model. At first, setup the layers:
    model <- keras_model_sequential()
    model %>%
      layer_dense(units = 64, activation = 'relu', input_shape = c(count_input_features)) %>%
      layer_dense(units = count.viruses_order, activation = 'softmax')

    # Compile the model:
    model %>% compile(
      optimizer = 'adam',
      loss = 'sparse_categorical_crossentropy',
      metrics = c('accuracy')
    )

    # Train the model:
    history <- model %>% fit(
      x_training, y_training, validation_data = list(x_validation, y_validation), verbose = 0,
      shuffle = T,
      epochs = 100,
      batch_size = floor(sqrt(nrow(x_training)))^2,
      callbacks = callback_early_stopping(monitor = "val_loss", patience = 10, verbose = 0, restore_best_weights = T)
    )

    # Evaluate accuracy:
    score <- model %>% evaluate(x_test, y_test)
    cat('Test loss:', score$loss, "\n")
    cat('Test accuracy:', score$acc, "\n")

    # Predict classes:
    predictions = model %>% predict_classes(x_test)
    (predictions = factor(predictions, levels = 0:(count.viruses_order-1)))

    Confusion.Matrix = data.frame(matrix(nrow = count.viruses_order, ncol = count.viruses_order+6))
    colnames(Confusion.Matrix) = c("TPR","TNR","PPV","NPV","Order","Order_ID",paste("Order",0:(count.viruses_order-1),sep = "_"))
    Confusion.Matrix$Order = names(count.viruses_all)
    Confusion.Matrix$Order_ID = 0:(count.viruses_order-1)
    for ( i in 0:(count.viruses_order-1) ) {
      my.indices = which(y_test == i )
      my.row = table(predictions[my.indices])
      Confusion.Matrix [i+1,-(1:6)] = my.row
    }
    total_sum = sum(Confusion.Matrix[,-(1:6)])
    for ( i in 0:(count.viruses_order-1) ) {
      my.matrix = Confusion.Matrix[,-(1:6)]
      TP = my.matrix [i+1,i+1]
      TN = sum(my.matrix [-(i+1),-(i+1)])
      C_P = sum(my.matrix[i+1,])
      C_N = sum(my.matrix[-(i+1),])
      P_C_P = sum(my.matrix[,i+1])
      P_C_N = sum(my.matrix[,-(i+1)])
      TPR = TP / C_P
      TNR = TN / C_N
      PPV = TP / P_C_P
      NPV = TN / P_C_N
      Confusion.Matrix [i+1,1:4] = round(100*c(TPR,TNR,PPV,NPV),1)

    }
    Confusion.Matrix

    accuracy_results[iteration] = score$acc
    confusion_matrix_results[[iteration]] = Confusion.Matrix
    acc_loss_graphics_results[[iteration]] = history

  }

  saveRDS(accuracy_results, paste0(temp.directory,"/accuracy_results.rds"))
  saveRDS(confusion_matrix_results, paste0(temp.directory,"/confusion_matrix_results.rds"))
  saveRDS(acc_loss_graphics_results, paste0(temp.directory,"/acc_loss_graphics_results.rds"))

}

#' Artificial neural network (ANN) classification of a new sample
#'
#' Step 2 of the taxa classification and estimation procedure.
#' This function loads an ANN model which was trained on artificial data,
#' computes input features of the new, adjusted sample file (FASTQ or SAM)
#' and stores the predicted classes of its read sequences.
#'
#' @seealso \code{\link{class_specific_PPVs}}, \code{\link{estimation_taxa_distribution}}
#' @import Biostrings
#' @import keras
#' @import stringi
#' @param fasta.file_path FASTA file path (e.g. FASTA file viral.genomic.fna downloaded from NCBI: ftp://ftp.ncbi.nih.gov/refseq/release/viral).
#' @param taxonomy.file_path Taxonomy file path (e.g. taxonomy file taxonomy_viruses_available.csv delivered within this package).
#' @param temp.directory Results directory containing accuracy results, generalised confusion matrix results and accuracy as well as loss graphics of the simulation runs.
#' @param read_sequences.file_path Only read sequences without identifier or species names, extracted from FASTQ or SAM file!
#' @param model.file_path File path of ANN model trained on artificially generated data.
#' @param predictions.file_path File path of the result file of predicted taxonomic orders.
#' @examples
#' # Please specify your file paths and directories!
#' fasta.file_path = "~/ag_bioinf/genomes/viruses_na/refseq/viral.genomic.fna" # Download from NCBI!
#' taxonomy.file_path = "inst/extdata/taxonomy_viruses_available.csv" # Relative file path
#' temp.directory = "~/ag_bioinf/research/metagenomics/temp" # Results directory
#' read_sequences.file_path = "~/ag_bioinf/research/metagenomics/Data/Seehund_Mapping/read_sequences.txt"
#' model.file_path = "inst/extdata/model_training_1_dataset.h5"
#' predictions.file_path = "~/ag_bioinf/research/metagenomics/temp/900000_training_samples_1_iteration/test/predictions_seal_sample.rds"
#' \dontrun{
#' ANN_new_sample ( fasta.file_path, taxonomy.file_path, temp.directory, read_sequences.file_path )
#' }
ANN_new_sample = function ( ) {
  read_sequences = readLines(read_sequences.file_path)
  # read_sequences = read_sequences [1:100]
  read_sequences.count = length(read_sequences)
  seq_len = unname(sapply(read_sequences,nchar))

  sequences = DNAStringSet(read_sequences)

  # One-mer, two-mer and three-mer distributions and inter-nucleotide distances (4 + 16 + 64 + 36 = 120 variables).
  two_mer_permutations = expand.grid(c("A","C","G","T"),c("A","C","G","T"))
  two_mer_permutations = as.matrix(two_mer_permutations)
  two_mer_permutations = apply(two_mer_permutations, MARGIN = 1, FUN = function (x) paste(x,collapse=""))
  two_mer_permutations = sort(two_mer_permutations)
  three_mer_permutations = expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T"))
  three_mer_permutations = as.matrix(three_mer_permutations)
  three_mer_permutations = apply(three_mer_permutations, MARGIN = 1, FUN = function (x) paste(x,collapse=""))
  three_mer_permutations = sort(three_mer_permutations)

  count_input_features = 4 + 16 + 64 + 36

  x_test = matrix(nrow = read_sequences.count, ncol = count_input_features)
  for ( i in 1:read_sequences.count) {
    if ( i %% 100 == 0 ) print(round(i/read_sequences.count*100,2))
    x_test [i,1:4] = letterFrequency(sequences[[i]], letters = c("A","C","G","T"), as.prob = F)
    for ( j in 5:20 ) {
      x_test [i,j] = stri_count_fixed(sequences[[i]],pattern = two_mer_permutations[j-4])
    }
    for ( j in 21:84 ) {
      x_test [i,j] = stri_count_fixed(sequences[[i]],pattern = three_mer_permutations[j-20])
    }
    dist_A = diff(stri_locate_all(sequences[[i]], regex = "A") [[1]] [,1])
    dist_A = table(factor(dist_A, levels = 2:10))
    dist_C = diff(stri_locate_all(sequences[[i]], regex = "C") [[1]] [,1])
    dist_C = table(factor(dist_C, levels = 2:10))
    dist_G = diff(stri_locate_all(sequences[[i]], regex = "G") [[1]] [,1])
    dist_G = table(factor(dist_G, levels = 2:10))
    dist_T = diff(stri_locate_all(sequences[[i]], regex = "T") [[1]] [,1])
    dist_T = table(factor(dist_T, levels = 2:10))
    x_test [i,85:93] = dist_A
    x_test [i,94:102] = dist_C
    x_test [i,103:111] = dist_G
    x_test [i,112:120] = dist_T
  }
  for ( my.index in 85:count_input_features ) {
    x_test[is.na(x_test[,my.index]),my.index] = 0
  }
  x_test = apply(x_test, MARGIN = 2, FUN = function(x) (x-min(x)) / diff(range(x)))
  colnames(x_test) = c("A","C","G","T",two_mer_permutations,three_mer_permutations,
                       paste0("d_A_",2:10),paste0("d_C_",2:10),paste0("d_G_",2:10),paste0("d_T_",2:10))

  # Load your previously on artificially generated reads trained ANN model ("e.g. model_training_1_dataset.h5"):
  model = load_model_hdf5(model.file_path)
  predictions = model %>% predict_classes(x_test)
  # Save the classification results of your model:
  saveRDS(predictions,predictions.file_path)
}

#' Prior and posterior estimation of taxa distribution
#'
#' Step 3 of the taxa classification and estimation procedure.
#' This function loads the predicted classes of the ANN model which was trained on artificial data,
#' computes prior as well as posterior taxa distribution estimations
#' and saves a graphics result file containing the estimation of the predicted classes.
#'
#' @seealso \code{\link{class_specific_PPVs}}, \code{\link{ANN_new_sample}}
#' @import Biostrings
#' @import keras
#' @import stringi
#' @import ggplot2
#' @import gridExtra
#' @import ggpubr
#' @import ggplotify
#' @param temp.directory Results directory containing accuracy results, generalised confusion matrix results and accuracy as well as loss graphics of the simulation runs.
#' @param graphics.directory Directory of graphics result file.
#' @param a_priori_table.file_path Predictions of previously classified reads of new sample file.
#' @examples
#' # Please specify your file paths and directories!
#' temp.directory = "inst/extdata/class_specific_PPVs_results/" # Results directory
#' graphics.directory = "~/ag_bioinf/research/metagenomics/ManuscriptNeuralNet/Graphics/"
#' a_priori_table.file_path = "inst/extdata/predictions_3154562.rds"
#' \dontrun{
#' estimation_taxa_distribution ( temp.directory, graphics.directory, a_priori_table.file_path )
#' }
estimation_taxa_distribution = function ( ) {
  estimation.a_priori.table = readRDS(a_priori_table.file_path)
  estimation.a_priori.table = estimation.a_priori.table[-10]

  confusion_matrix_results = readRDS(paste0(temp.directory,"confusion_matrix_results.rds"))
  loss_acc_results = readRDS(paste0(temp.directory,"loss_acc_results.rds"))
  loss_acc_graphics_results = readRDS(paste0(temp.directory,"loss_acc_graphics_results.rds"))

  order_names = c("Bunyavirales","Caudovirales","Herpesvirales","Ligamenvirales","Mononegavirales",
                  "Nidovirales","Ortervirales","Picornavirales","Tymovirales")
  TPR = matrix(NA_real_, nrow = 9, ncol = 10)
  PPV = matrix(NA_real_, nrow = 9, ncol = 10)
  for ( i in 1:10 ) {
    TPR[,i] = confusion_matrix_results [[i]] [,1]
    PPV[,i] = confusion_matrix_results [[i]] [,3]
  }
  rownames(TPR) = rownames(PPV) = order_names
  df.TPR_and_PPV = data.frame(Order = rep(order_names, each = 10), TPR = as.vector(t(TPR)), PPV = as.vector(t(PPV)))

  estimation.a_priori = as.data.frame(matrix(nrow = 9, ncol = 12))
  colnames(estimation.a_priori) = c("Order","Order_ID,",paste0("Iteration_",1:10))
  estimation.a_priori[,1] = order_names
  estimation.a_priori[,2] = 0:8
  for ( j in 3:12 ) {
    estimation.a_priori[,j] = unlist(estimation.a_priori.table)
  }

  estimation.a_posteriori = estimation.a_priori


  PPV_probabilities = vector("list",10)
  for ( j in 1:10 ) {
    PPV_probabilities[[j]] = matrix(NA_real_,9,9)
    rownames(PPV_probabilities[[j]]) = paste0("Order_",0:8)
    colnames(PPV_probabilities[[j]]) = paste0("Order_",0:8)
  }
  for ( i in 1:9 ) {
    for ( j in 1:10 ) {
      PPV_probabilities[[j]] [,i] = confusion_matrix_results[[j]] [,i+6] / sum(confusion_matrix_results[[j]] [,i+6])
    }
  }

  PPV_probabilities.mean = vector("list",10)
  for ( j in 1:10 ) {
    PPV_probabilities.mean[[j]] = matrix(NA_real_,9,9)
    rownames(PPV_probabilities.mean[[j]]) = paste0("Order_",0:8)
    colnames(PPV_probabilities.mean[[j]]) = paste0("Order_",0:8)
    for ( i in 1:9 ) {
      temp = sapply(PPV_probabilities[-j], FUN = function(x) x[,i])
      PPV_probabilities.mean[[j]] [,i] = rowMeans(temp)
    }
  }

  for ( j in 1:10 ) {
    estimation.a_posteriori[,j+2] = PPV_probabilities.mean[[j]] %*% matrix(estimation.a_priori[,j+2],ncol = 1)
    estimation.a_posteriori[,j+2] = PPV_probabilities[[j]] %*% matrix(estimation.a_priori[,j+2],ncol = 1)
  }

  my.data_frame = data.frame(Order = rep(rep(order_names, each = 10),2), Estimation = c(rep("Prior",9*10),rep("Posterior",9*10)),
                             Iteration = rep(1:10, 2*9), Count = c(as.vector(t(as.matrix(estimation.a_priori[,-(1:2)]))),
                                                                   as.vector(t(as.matrix(estimation.a_posteriori[,-(1:2)])))))
  my.data_frame$Order = as.factor(my.data_frame$Order)
  my.data_frame$Estimation = factor(my.data_frame$Estimation, levels = c("Prior","Posterior"))
  my.data_frame$Iteration = as.factor(my.data_frame$Iteration)

  my.data_frame.relative = my.data_frame
  for ( i in 1:10 ) {
    my.data_frame.relative[seq(0+i,90,by=10),4] = my.data_frame.relative[seq(0+i,90,by=10),4] / sum(my.data_frame.relative[seq(0+i,90,by=10),4]) * 100
    my.data_frame.relative[seq(90+i,180,by=10),4] = my.data_frame.relative[seq(90+i,180,by=10),4] / sum(my.data_frame.relative[seq(90+i,180,by=10),4]) * 100
  }
  my.data_frame.relative$Order = factor(my.data_frame.relative$Order, levels = rev(c("Bunyavirales","Caudovirales","Herpesvirales","Ligamenvirales","Mononegavirales",
                                                                                     "Nidovirales","Ortervirales","Picornavirales","Tymovirales") ))

  require("ggplot2")
  .df <- data.frame(x = my.data_frame.relative$Order, y =
                      my.data_frame.relative$Count, z = my.data_frame.relative$Estimation)
  .plot <- ggplot(data = .df, aes(x = factor(x), y = y, colour = z)) +
    stat_summary(fun.y = "mean", geom = "point",  position =
                   position_dodge(width = 0.6) ) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",  position =
                   position_dodge(width = 0.6), pch = 10, size = 1,  width = 0.1, fun.args = list(conf.int =
                                                                                                    1.0)) +
    coord_flip() +
    xlab("Order") +
    ylab("Relative frequency (%)") +
    labs(colour = "Estimation") +
    theme_bw(base_size = 25, base_family = "sans")
  print(.plot)

  ggsave(filename =
           paste0(graphics.directory,"prior_posterior_estimation.pdf"),
         plot = .plot, width = 12, height = 8)

  rm(.df, .plot)
  dev.off()

}
