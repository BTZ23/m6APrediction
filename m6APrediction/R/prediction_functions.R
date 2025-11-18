#' DNA Sequence Encoding for m6A Prediction
#'
#' This function converts DNA sequences into a feature matrix for machine learning.
#'
#' @param dna_strings A character vector of DNA sequences (e.g., c("GGACA", "ACGUU"))
#' @return A data.frame with encoded nucleotide positions as factors
#' @examples
#' dna_encoding(c("GGACA", "ACGUU"))
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Batch m6A Site Prediction
#'
#' Predict m6A modification sites for multiple sequences using a trained model.
#'
#' @param ml_fit A trained machine learning model (e.g., randomForest)
#' @param feature_df A data.frame containing features for prediction. Must include columns:
#'   gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, 
#'   evolutionary_conservation, DNA_5mer
#' @param positive_threshold Probability threshold for positive classification (default: 0.5)
#' @return A data.frame with original features plus predicted probability and status
#' @examples
#' \dontrun{
#' # Load model and data first
#' predictions <- prediction_multiple(rf_model, my_features)
#' }
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", 
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  
  dna_encoded <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df, dna_encoded)
  
  predictions <- predict(ml_fit, newdata = feature_df, type = "prob")
  feature_df$predicted_m6A_prob <- round(predictions[, "Positive"], 3)
  feature_df$predicted_m6A_status <- ifelse(feature_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")
  
  return(feature_df)
}

#' Single m6A Site Prediction
#'
#' Predict m6A modification for a single sequence.
#'
#' @param ml_fit A trained machine learning model
#' @param gc_content Numeric GC content value (0-1)
#' @param RNA_type Character RNA type ("mRNA", "lincRNA", "lncRNA", "pseudogene")
#' @param RNA_region Character RNA region ("CDS", "intron", "3'UTR", "5'UTR")
#' @param exon_length Numeric exon length
#' @param distance_to_junction Numeric distance to junction
#' @param evolutionary_conservation Numeric conservation score (0-1)
#' @param DNA_5mer Character 5-mer DNA sequence (e.g., "GGACA")
#' @param positive_threshold Probability threshold for positive classification (default: 0.5)
#' @return A named vector with predicted probability and status
#' @examples
#' \dontrun{
#' result <- prediction_single(rf_model, 0.5, "mRNA", "CDS", 10, 8, 0.5, "GGACA")
#' }
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, 
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  
  result_df <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  
  returned_vector <- c(
    predicted_m6A_prob = as.character(result_df$predicted_m6A_prob),
    predicted_m6A_status = result_df$predicted_m6A_status
  )
  
  return(returned_vector)
}

