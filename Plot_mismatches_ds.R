#!/usr/bin/env Rscript
# 2021-03-26
# Becky Cribdon and Rosie Everett
#
# Makes MetaDamage plots from ".mismatches.txt" files.
# Also outputs 95% confidence intervals for C->T and G->A mismatches at the end positions as a .txt table.

library(ggplot2) # Plots
library(gridExtra) # grid.arrange
library(grid) # textGrob

output_directory <- "MetaDamage_outputs/"

line_colour_CT <- "red" # C->T.
line_colour_GA <- "blue" # G->A.
line_colour_other <- "grey" # All other substitutions. C-G, C-A, T-C, T-G, T-A, G-C, G-T, A-C, A-T, A-G.

mismatch_files <- commandArgs(trailingOnly=TRUE)
#mismatch_files[1] <- "MetaDamage_outputs/test1.mismatches.txt" # TESTING ONLY

#number_of_samples <- length(mismatch_files) # Only one input file is expected.

# Generate the output file names:
mismatch_basenames <- basename(mismatch_files) # Remove parent directories.
mismatch_basenames <- substr(mismatch_basenames[1], 1, (nchar(mismatch_basenames[1])-15) ) # Remove '.mismatches.txt'.
output_name_plot <- paste(output_directory, mismatch_basenames[1], '.MetaDamage_plots.pdf', sep='') # A PDF for the plot.
output_name_CIs <- paste(output_directory, mismatch_basenames[1], '.MetaDamage_CIs.txt', sep='') # A .txt for the CI information.

cat("Input\tC->T at 5' end\tC->T 95% CI lower\tC->T 95% CI upper\tG->A at 3' end\tG->A 95% CI lower\tG->A 95% CI upper\n", file=output_name_CIs) # Write a header to the output text file.


# Generate left plot
#-------------------
make_left_plot <- function(i) { # Make the plotting via functions so we can do it over and over.
  count_stats <- readLines(con=mismatch_files[i], n=2) # Read in the top two lines.
  
  n_reads <- count_stats[1]
  n_reads <- as.numeric( substr(n_reads, 30, nchar(n_reads)) ) # Remove the text label and treat the remaining digit(s) as a number.
  
  n_trials <- count_stats[2] # The number of Cs at position 0 in the reference sequences.
  n_trials <- as.numeric( substr(n_trials, 34, nchar(n_trials)) ) # Remove the text label and treat the remaining digit(s) as a number.
  
  mismatches = read.table(mismatch_files[i], skip=4, nrows=25, sep="\t", stringsAsFactors=F, quote="\"") # Just take the first half of the table, excluding header.
  colnames(mismatches) <- c("Position","P_AtoT","P_AtoC","P_AtoG","P_TtoA","P_TtoC","P_TtoG","P_CtoA","P_CtoT","P_CtoG","P_GtoA","P_GtoT","P_GtoC") # R-friendly characters.
  
  if (n_reads == 0) {
    position_0_CT_mismatch <- NA
    CI_lower_bound <- NA
    CI_upper_bound <- NA
    is_over_0.3 <- 0
  } else {
    # Calculate 95% CI bounds for C->T at position 0:
    position_0_CT_mismatch <- mismatches$P_CtoT[1]
    alpha <- (n_trials * position_0_CT_mismatch) + 1 # Number of successes + 1.
    beta <- (n_trials * (1-position_0_CT_mismatch)) + 1 # Number of failures + 1.
    CI_lower_bound = 0 # Default the lower bound to 0; it doesn't make sense for it to go below 0.
    if (position_0_CT_mismatch > 0) {
      CI_lower_bound <- qbeta(0.025, shape1=alpha, shape2=beta) # pbeta() is the cumulative distribution function. qbeta() is like the inverse, calculating the x (mismatch value) given a y (probability).
    }
    CI_upper_bound <- qbeta(0.975, shape1=alpha, shape2=beta)
    
    is_over_0.3 <- 0 # A flag.
    if ( (!is.na(CI_upper_bound)) && (CI_upper_bound > 0.3) ) { # First check the upper CI bound. If it's >0.3,
      is_over_0.3 <- 1 # Set the flag to 1.
    } else if (!(is.na(mismatches$P_AtoT[1]))) { # Then check the mismatch values. If the first value in the table isn't NA (which means all are):
      for (column in (2:ncol(mismatches)) ) { # Check every data column for values >0.3.
        if (! all(mismatches[,column] < 0.3)) { # If any mismatches value is >0.3, 
          is_over_0.3 <- 1 # Set the flag to 1.
        }
      }
    }
  }
  
  cat(paste(mismatch_basenames[i], position_0_CT_mismatch, CI_lower_bound, CI_upper_bound, sep="\t"), file=output_name_CIs, append=T) # Add the input name and CI bounds to the output file.
  
  baseplot <- ggplot(mismatches) + # Set up a base plot with most of the formatting.
    xlab("Position from 5' end") +
    ylab(expression('P'[substitution])) +
    scale_x_continuous(breaks=c(0:24),
                       limits=c(-0.25,24),
                       labels=factor(0:24), # Make it label every base in use.
                       expand=c(0,0.5) # Reduce padding between the data and axis.
                       ) + 
    scale_y_continuous( breaks=c(seq(from=0, to=0.3, by=0.05)),
          limits=c(0, 0.3),
          labels=c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"), # Specify the breaks to 2dp to maintain the same y-axis width as plots that go up to 1.
          expand=c(0.025,0) # Reduce padding between the data and axis.
          ) +
    theme_bw() +
    theme(aspect.ratio=4/6.5,
          axis.title=element_text(size=6),
          axis.title.x=element_text(margin=margin(t=6.25, r=0, b=0, l=0)), # Hack to account for lack of '-' compared to right plot.
          axis.text=element_text(size=6),
          axis.text.x=element_text(angle=90),
          axis.line=element_line(colour="black"),
          panel.border=element_blank(),
          panel.grid=element_blank(), # Remove gridlines.
          plot.margin = unit(c(0,5.5,5.5,5.5), "pt") # Remove the top margin. Leave the rest on default.
          )

  # Now modify the base plot depending on the data:
  if (n_reads == 0) { # If there weren't any reads, fill the plot with "No data".
      plot <- baseplot +
        geom_text(x=12.5, y=0.15, size=3, label="No data") 
  } else {
    if (is_over_0.3 == 1) { # If the maximum value was over 0.3:
        cat(paste("\t", mismatch_basenames[i], ": high 5'-3' mismatch values. Increasing y-axis limit to 1.\n", sep=""))
        suppressMessages( # Wrap this in suppressMessages() to omit "Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.".
          baseplot <- baseplot +
          geom_hline(yintercept = 0.3, size=0.5, color="grey", linetype="dashed") + # Add a further line at 0.3 to help compare with other plots.
          ylim(0, 1) # Make axis limits constant for all particularly-high plots.
        )
    }
    plot <- baseplot + # Plot the data with the two most important on top:
        geom_line(aes(x=Position, y=P_AtoT), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_AtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_AtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_CtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_CtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoT), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoA), size=0.5, col=line_colour_GA) +
        geom_line(aes(x=Position, y=P_CtoT), size=0.5, col=line_colour_CT)
    if (!is.na(CI_lower_bound)) { # If there was a CI, plot it as error bars.
      plot <- plot + geom_errorbar(aes(x=0, ymin=CI_lower_bound, ymax=CI_upper_bound), width=0.5, size=0.25)
    }
    plot <- plot + geom_point(x=0, y=position_0_CT_mismatch, size=0.5, col=line_colour_CT)
    plot # Print the plot.
  }
} # End of function.


make_right_plot <- function(i) {
  count_stats <- readLines(con=mismatch_files[i], n=3) # Read in the top three lines.
  n_reads <- count_stats[1]
  n_reads <- as.numeric( substr(n_reads, 30, nchar(n_reads)) ) # Remove the text label and treat the remaining digit(s) as a number.

  n_trials <- count_stats[3] # The number of Gs at position 0 in the reference sequences.
  n_trials <- as.numeric( substr(n_trials, 34, nchar(n_trials)) )
  
  mismatches = read.table(mismatch_files[i], skip=29, header=F, sep="\t", stringsAsFactors=F, quote="\"") # Just take the second half of the table (therefore no header).
    colnames(mismatches) <- c("Position","P_AtoT","P_AtoC","P_AtoG","P_TtoA","P_TtoC","P_TtoG","P_CtoA","P_CtoT","P_CtoG","P_GtoA","P_GtoT","P_GtoC") # R-friendly characters.
  
  if (n_reads == 0) {
    position_0_GA_mismatch <- NA
    CI_lower_bound <- NA
    CI_upper_bound <- NA
    is_over_0.3 <- 0
  } else {
    # Calculate 95% CI bounds for G->A at position 0:
    position_0_GA_mismatch <- mismatches$P_GtoA[25]
    alpha <- (n_trials * position_0_GA_mismatch) + 1 # Number of successes + 1.
    beta <- (n_trials * (1-position_0_GA_mismatch)) + 1 # Number of failures + 1.
    CI_lower_bound <- 0 # Default the lower bound to 0; it doesn't make sense for it to go below 0.
    if (position_0_GA_mismatch > 0) {
      CI_lower_bound <- qbeta(0.025, shape1=alpha, shape2=beta) # pbeta() is the cumulative distribution function. qbeta() is like the inverse, calculating the x (mismatch value) given a y (probability).
    }
    CI_upper_bound <- qbeta(0.975, shape1=alpha, shape2=beta)
    
    is_over_0.3 <- 0 # A flag.
    if ( (!is.na(CI_upper_bound)) && (CI_upper_bound > 0.3) ) { # First check the upper CI bound. If it's >0.3,
      is_over_0.3 <- 1 # Set the flag to 1.
    } else if (!(is.na(mismatches$P_AtoT[1]))) { # Then check the mismatch values. If the first value in the table isn't NA (which means all are):
      for (column in (2:ncol(mismatches)) ) { # Check every data column for values >0.3.
        if (! all(mismatches[,column] < 0.3)) { # If any mismatches value is >0.3, 
          is_over_0.3 <- 1 # Set the flag to 1.
        }
      }
    }
  }
  
  cat(paste("\t", position_0_GA_mismatch, "\t", CI_lower_bound, "\t", CI_upper_bound, "\n", sep=""), file=output_name_CIs, append=T) # Add the CI bounds to the output file.
  
  baseplot <- ggplot(mismatches) + # Set up a base plot with most of the formatting.
    xlab("Position from 3' end") +
    ylab(expression('P'[substitution])) +
    #xlim(-24,0)
    scale_x_continuous( breaks=c(-24:0),
                       limits=c(-24,0.5),
                       labels=factor(-24:0), # Make it label every base in use.
                       expand=c(0,0.5) # Reduce padding between the data and axis.
                       ) + 
    scale_y_continuous( breaks=c(seq(from=0, to=0.3, by=0.05)),
          limits=c(0, 0.3),
          labels=c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"), # Specify the breaks to 2dp to maintain the same y-axis width as plots that go up to 1.
          expand=c(0.025,0), # Reduce padding between the data and axis.
          position="right"
          ) +
    theme_bw() +
    theme(aspect.ratio=4/6.5,
          axis.title=element_text(size=6),
          axis.text=element_text(size=6),
          axis.text.x=element_text(angle=90),
          axis.line=element_line(colour="black"),
          panel.border=element_blank(),
          panel.grid=element_blank(), # Remove gridlines.
          plot.margin = unit(c(0,5.5,5.5,5.5), "pt") # Remove the top margin. Leave the rest on default.
          )

  # Now modify the base plot depending on the data:
  if (n_reads == 0) { # If there weren't any reads, fill the plot with "No data".
      plot <- baseplot +
        geom_text(x=-12.5, y=0.15, size=3, label="No data") 
  } else {
    if (is_over_0.3 == 1) { # If the maximum value was over 0.3:
        cat(paste("\t", mismatch_basenames[i], ": high 3'-5' mismatch values. Increasing y-axis limit to 1.\n", sep=""))
        suppressMessages( # Wrap this in suppressMessages() to omit "Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.".
          baseplot <- baseplot +
          geom_hline(yintercept = 0.3, size=0.5, color="grey", linetype="dashed") + # Add a further line at 0.3 to help compare with other plots.
          scale_y_continuous( breaks=c(seq(from=0, to=1, by=0.1)),
          limits=c(0, 1),
          labels=c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90", "1.00"), # Specify the breaks to 2dp to maintain the same y-axis width as plots that go up to 1.
          #expand=c(0.025,0), # Reduce padding between the data and axis.
          position="right"
          ) # Make axis limits constant for all particularly-high plots.
        )
    }
    plot <- baseplot + # Plot the data with the two most important on top:
        geom_line(aes(x=Position, y=P_AtoT), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_AtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_AtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_CtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_CtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_TtoG), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoA), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoT), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoC), size=0.5, col=line_colour_other) +
        geom_line(aes(x=Position, y=P_GtoA), size=0.5, col=line_colour_GA) +
        geom_line(aes(x=Position, y=P_CtoT), size=0.5, col=line_colour_CT)
    if (!is.na(CI_lower_bound)) { # If there was a CI, plot it as error bars.
      plot <- plot + geom_errorbar(aes(x=0, ymin=CI_lower_bound, ymax=CI_upper_bound), width=0.5, size=0.25)
    }
    plot <- plot + geom_point(x=0, y=position_0_GA_mismatch, size=0.5, col=line_colour_GA)
      
    if (is_over_0.3 == 1) { # Add the read count to the top right of the plot. hjust=1 justifies it right, so it should always be inside the box.
      plot <- plot + geom_text(x=-0.5, y=0.99, size=2, hjust=1, colour="grey30", label=n_reads)
    } else {
      plot <- plot + geom_text(x=-0.5, y=0.29, size=2, hjust=1, colour="grey30", label=n_reads)
    }
    plot # Print the plot.
  }
} # End of function.

list_of_left_and_right_plots <- list() # Temporarily holds a pair of plots.
list_of_paired_plots <- list() # Holds all pairs.

#for (i in (1:number_of_samples)) {
  i=1
  list_of_left_and_right_plots[1] <- lapply(i, make_left_plot)
  list_of_left_and_right_plots[2] <- lapply(i, make_right_plot)
  paired_plot <- grid.arrange(
    list_of_left_and_right_plots[[1]],
    list_of_left_and_right_plots[[2]],
    ncol=2,
    top=textGrob(mismatch_basenames[i], gp=gpar(fontsize=8), x=0, hjust=0) # The hjust justifies it left, so it won't go off the edge.
    )
  list_of_paired_plots[[i]] <- paired_plot
#}

#pdf(file=output_name_plot, width=6, height=2.2*number_of_samples)
pdf(file=output_name_plot, width=6, height=2.2)
grid.arrange(grobs=list_of_paired_plots, ncol=1) # Print the plots to it.
invisible(dev.off()) # Close the PDF.

