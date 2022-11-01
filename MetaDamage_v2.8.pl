#!/usr/bin/perl
use strict;
use warnings;

# Rosie Everett, Becky Cribdon

#To run: perl MetaDamage_v2.8.pl -fastas [input_fasta] -threads [number of processing threads] 

############### set up ############################

# Set up of flags 
use Getopt::Long;

# Variables defined 
my $fastas = '';
my $threads = '';
my $blast_db = ''; # Default the local BLAST database to empty, needs to be defined 
my @fas = '';
my $fasta = '';
my $fasta_filename = '';
my $batches = 100;
my $batch = 0;
my $batchcount = 0;
my $batchround = 0;
my @reference = '';
my $reference_ID = (); 
my $reference_sequence = ();
my $entry = ();
my $reference = '';
my $qseqid = '';
my @fstarts = ();
my @fstops = ();
my @fqseqids = ();
my @bigfqseqids = ();
my @bigfstarts = ();
my @bigfstops = (); 
my @fsseqids = ();
my @bigfsseqids = ();
my @freversed = ();
my @bigfreversed = ();
my @fquery_seqs = ();
my @bigfquery_seqs = ();
my $fsseqids = ();
my $bigcheck = 0;
my $noblast = 0;
my $blastfile = '';
my $basename = '';
my $single_strand = '';
my $single_s = 0;
my $help = '';
my $outfile = '';
my @blasts = '';
my $counter = 0;
GetOptions (
	'fastas=s' => \$fastas,	#comma-separated list of input FASTA files (.fasta, .fna etc)
	'threads=i' => \$threads, #number of threads used to process analysis 
	'blast_db=s' => \$blast_db, # If the option is called, overwrite the default with the option input.
	'batches=i' => \$batches, #set up batch sizes to processes sequence groups
	'blastfile=s' => \$blastfile, #option to use a blast output prepared already
	'single_strand=s' => \$single_strand, #option to return damage profile for single stranded library prep recording C to T at both ends value TRUE  ss
	'help=s' => \$help, # option to display flag information
	'out=s' => \$outfile, # option to change the basename of output files 
);

if ($help) {
	printf "\t################### MetaDamage v2.8 ######################
\t USAGE: perl MetaDamage_v2.8.pl -flag [options]
\t FLAGS
\t -fastas : fasta input file, essential for all options
\t -threads : number of threads for blast algorithm to use
\t -blast_db : database to search. Optional, will default to searching NCBI if no database specified. This is SLOW!
\t -batches : batch size of sequences to run in analysis to make processing more efficient for high numbers of queries, default batch size 100
\t -blastfile : option to input a prepared blast output and skip blast search
\t -single_strand : TRUE to generate mismatch profiles suitable for interpretation with single stranded library preps
\t -out : option for output files root name
\t -help : help option to generate this output.
";
exit;

}
if ($single_strand) {
	$single_s = 1;
	printf "\n\n\tCalculating for single stranded libraries....\n";
}
if ($blastfile) {
		$noblast = 1;
		print "\n\n\tWorking from a ready blast output file.\n\tNote you must have used format flag -outfmt \"6 std qlen\"\n";
		@blasts = split /\,/, $blastfile;
		unless ($fastas) { die "\nA fasta input file is still required,\nspecify a comma-separated list of input FASTA files (.fasta, .fna, etc.) using -fastas\n";}
}

# Check the necesssary variables were given before we do anything else:
unless ($fastas) { # Note: you could use a regex for this, but regexes are slow, so use alternatives where possible.
    die "\nSpecify a comma-separated list of input FASTA files (.fasta, .fna, etc.) using -fastas\n";
}
if ($threads eq '' or $threads == 0) { # Zero threads is not valid either.
    unless ($blastfile) { die "\nSpecify a number of threads using -threads\n";}
}
unless ($blast_db) {
    #die "\nSpecify a local BLAST database using -blast_db\n";
	print "\n\n\n\tNo local database specified, defaulting to remote settings and NCBI E-utilities. Warning: This will take longer!\n\n\n";

}

my $error_log_filename = 'efetch_error_log.txt';

printf "\t################ Welcome to MetaDamage analysis v2.8 ################\n";
printf "\tCreating output directory\n";
my $output_directory = 'MetaDamage_outputs';
mkdir $output_directory;  # Make a master output folder to collect all output in.

############### Step 1: BLAST ############################
### Run BLAST 
@fas = split /\,/, $fastas;
foreach $fasta (@fas) {
    
    my @fasta = split('\.', $fasta); # Backslash to escape the special properties of '.'.
if ($outfile eq '') {   $basename = join('.', @fasta[0..($#fasta-1)]); # Join with a '.' elements 0 to (end-1) of @fasta.
			} else {$basename = $outfile;}
    my $output = $output_directory . '/' .$basename .".blast.txt";
    if ($noblast == 0) {
    printf "\tRunning BLAST to identify reference sequences...\n";

    if ($blast_db ne '') {
   `blastn -db $blast_db -num_threads $threads -query $fasta -out $output -max_target_seqs 1 -max_hsps 1 -outfmt \"6 std qlen\"`;	  
    } else {
	print "\tRemote BLAST search.....\n";
   `blastn -db nt -remote -query $fasta -out $output -max_target_seqs 1 -max_hsps 1 -outfmt \"6 std qlen\"`;
    }
    } else { print "\tUsing input blast file....\n";}

    ############### Step 2: Retrieve reference files #############
    my $blast_filename = ();
if ($noblast == 1) {
	$blast_filename = $blasts[$counter];
} else {
       $blast_filename = $output_directory . '/' ."$basename.blast.txt";
}
    my $paired_filename = $output_directory . '/' .$basename . '.paired.txt';

    printf "\tRetrieving reference files...\n";
    open (my $paired_filehandle, '>', $paired_filename) or die "\t\tCould not open output file $paired_filename for writing: $!\n"; # Open the paired output file.
    



    # Store the whole BLAST file in a hash. The FASTA will probably be bigger, so we will work through that line by line instead of the other way around.
	my @qseqid = (); my @query_seqs = (); my $qseq = (); my @reversed = (); 

    if (! open(my $blast_filehandle, $blast_filename) ) { # If the BLAST file won't open successfully, assume it's missing and skip this FASTA.
        print "\t\tERROR: Could not open BLAST file $blast_filename: $!\n";
        exit;
    } else {
        my %blast_file = (); # Keys are (unique) query sequence IDs, values are the rest of the BLAST entry.

	if ($batches > 0) {
		$batch = 0; $batchcount = 1;
		open (BATCH, ">batchfile.$batchcount") or die "Cannot write batchfile\n";
	}
        while (1) { # Run this loop until "last" is called.
            my $line = readline($blast_filehandle);
            if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
            
            chomp $line;
            my @line = split("\t", $line);
            $blast_file{$line[0]} = join("\t", @line[1..$#line]);
        }
        close $blast_filehandle;
        
        if (!keys %blast_file) { exit; } # If the BLAST hash is empty, either because the BLAST file was empty or there wasn't one, finish the script here.
        
    
        # Open the FASTA and read record by record.
        open(my $fasta_filehandle, $fasta) or die "\t\tCould not open original FASTA $fasta_filename: $!\n";
        
        $/ = '>'; # Set the record separator to '>', which separates FASTA records.
        readline $fasta_filehandle; # Skip the first 'record'. That's just the first separator: a single '>'.

        my $count_reads_discarded = 0; # A count of reads discarded because they do not align with the start of their reference.
    	my $first_efetch_start = 1; # A flag to say whether efetch has been called before or not.

        
        while (1) { # Run this loop until "last" is called.
                
                my $record = readline($fasta_filehandle);
                if (! defined $record) { last }; # If there is no next record, exit the loop. You've processed the whole file.
                
                my @record = split("\n", $record); # Split the record into its two lines.
		$qseqid = '';
                $qseqid = $record[0]; # Line 0 is the ID.
                chomp $qseqid;
                
                
                # Search for the query name and in the BLAST hash.
                if (exists $blast_file{$qseqid}) {
                    
                    my $blast_entry = "$qseqid\t" . $blast_file{$qseqid}; # Key\tvalue.
                    my @blast_entry = split("\t", $blast_entry); # Split the lot into columns.
    
                    my $sstart = $blast_entry[8];
                    my $send = $blast_entry[9];
                    my $reversed = 0; # Default reference sequence direction to the correct way round (0).
                    my $new_start = 0;
                    my $new_end = 0;
                    
                    if ($sstart < $send) { # If the reference matches from 5'-3' (same direction as the query),
                        $new_start = $sstart - ($blast_entry[6]-1); # New start coordinate = [start coordinate of the alignment in the reference] - ([start coordinate of the alignment in the query]-1).
                        $new_end = $send + ($blast_entry[12]-$blast_entry[7]); # $new_end is [end coordinate of the alignment in the reference] + ([query length]-[end coordinate of the alignment in the query]).
                    } elsif ($sstart > $send) { # If the reference matches from 3'-5' (opposite direction to the query),
                        $reversed = 1; # Mark as reversed.
                        $new_start = $send - ($blast_entry[12]-$blast_entry[7]);
                        $new_end = $sstart + ($blast_entry[6]-1);
                    } else {
                        print "\t\tError in BLAST file?\n";
                        exit;
                    }
                    
                    if ($new_start < 1) { # The start of the reference sequences is always position 1. If you have a start <1, the start of the query does not align with the reference, so we can't get any signal from the 5' end. Exclude these.
                        $count_reads_discarded ++;
                        next;
                    }
    
                    # Look up the reference sequence using blastcmd and extract it according to the new coordinates.
                    my $sseqid = $blast_entry[1];
                    $reference = ''; # Because the efetch command is run in backticks, its output can be saved as an object. It is a line-wrapped FASTA followed by a newline.
        
                    my $range = "$new_start" . '-' . "$new_end"; # blastdbcmd requires the start and end coordinates in the format x-y.
		    if ($batches > 0) {
			$qseq = $record[1];
			push (@qseqid, $qseqid); push (@query_seqs, $qseq); push (@reversed, $reversed);
			print BATCH "$sseqid $new_start-$new_end\n";
			# and collect query info too here in batches
			unless ($blast_db ne '') {
				# prepare 'batches' for efetch
				#print "start coordinate: $new_start\n";
				if ($new_start < 100000) {
						# pick up non genome sized entries here for 'batching'
					push (@fstarts, $new_start);
					push (@fstops, $new_end);
					push (@fqseqids, $qseqid);		
					push (@fsseqids, $sseqid);
					push (@freversed, $reversed);
					push (@fquery_seqs, $qseq);
				} else {
						# very large (genome) sized entries to be dealt with individually
                                        push (@bigfstarts, $new_start);
                                        push (@bigfstops, $new_end);
                                        push (@bigfqseqids, $qseqid);
                                        push (@bigfsseqids, $sseqid);
                                        push (@bigfreversed, $reversed);
					push (@bigfquery_seqs, $qseq);
					$bigcheck = 1;
				}
			} 
			$batch++;
			if ($batch == $batches) {
			close BATCH; $batch = 0;
			# pack up and ship out
			if ($blast_db ne '') {
			$reference = `blastdbcmd -db $blast_db -entry_batch batchfile.$batchcount`;

			if ($reference eq '') {
                        print "\n\nERROR: blastdbcmd cannot retrieve sseqid batchfile.$batchcount\n"; exit;
                    	}
			} else {
				# efetch
				$fsseqids = join ("\,", @fsseqids);
                                        while (1) { # Run this loop until "last" is called.
#                                                $reference = `efetch -db Nucleotide -id $fsseqids -format fasta 2>> $error_log_filename`; # The efetch command.
                                                 #print "\n\nEFETCH COMMAND: efetch -db nuccore -id $fsseqids -format fasta 2>> $error_log_filename";
                                                 $reference = `efetch -db nuccore -id $fsseqids -format fasta 2>> $error_log_filename`; # The efetch command.


                                                if ($reference eq '') { # Failed efetches return a newline. If the return is completely empty, then efetch couldn't start for some reason - suggesting a different error.
                                                if ($first_efetch_start == 1) {
                                                print "\t\tERROR: efetch could not start. See $error_log_filename.\n";
                                                exit;
                                                } else {
                                                print "\t\tERROR: efetch could not start for sequence $qseqid. See $error_log_filename.\n";
                                                exit;
                                                }
                                        }

                                my @reference = split("\n", $reference);
                                if ( ($reference[1])  && ($reference =~ /\S/) ) { # If the reference was successfully split and contains at least one non-whitespace character, exit the loop.
                                $first_efetch_start = 0;
                                last;
                                }
                                }

			}
			unlink ("batchfile.$batchcount");
			$batchcount++;
			$batchround = 1;
			open (BATCH, ">batchfile.$batchcount") or die "Cannot write batchfile.$batchcount\n";
			}
			} else {
			    if ($blast_db ne '') {
                    		$reference = `blastdbcmd -db $blast_db -entry $sseqid -range $range`;
                       
                    		if ($reference eq '') {
                        		print "\n\nERROR: blastdbcmd cannot retrieve sseqid '$sseqid' (qseqid: '$qseqid')\n"; exit;
                    		}
			    } else {
				# default to efetch individual
                    			while (1) { # Run this loop until "last" is called.
#                        			$reference = `efetch -db Nucleotide -id $sseqid -seq_start $new_start -seq_stop $new_end -format fasta 2>> $error_log_filename`; # The efetch command.
                                                #print "\n\nEFETCH COMMAND: efetch -db nuccore -id $sseqid -seq_start $new_start -seq_stop $new_end -format fasta 2>> $error_log_filename";
                                                $reference = `efetch -db nuccore -id $sseqid -seq_start $new_start -seq_stop $new_end -format fasta 2>> $error_log_filename`; # The efetch command.
                        			if ($reference eq '') { # Failed efetches return a newline. If the return is completely empty, then efetch couldn't start for some reason - suggesting a different error.
                            			if ($first_efetch_start == 1) {
                                		print "\t\tERROR: efetch could not start. See $error_log_filename.\n";
                                		exit;
                            			} else {
                                		print "\t\tERROR: efetch could not start for sequence $qseqid. See $error_log_filename.\n";
                                		exit;
                            			}
                        		}

                        	my @reference = split("\n", $reference);
                        	if ( ($reference[1])  && ($reference =~ /\S/) ) { # If the reference was successfully split and contains at least one non-whitespace character, exit the loop.
                            	$first_efetch_start = 0;
                            	last;
                        	}
                    		}
				
			    }
		    }	
                unless ($batches > 0) {    
                    # Separate the header and put all sequence lines together.
                    @reference = split("\n", $reference);
                    $reference_ID = $reference[0];
                    $reference_sequence = join('', @reference[1..$#reference]);
    
                    # Is the reference section 3'-5'? If so, reverse complement it.
                    if ($reversed == 1) {
                        $reference_sequence = reverse($reference_sequence);
                        $reference_sequence =~ tr/[ACGTacgt]/[TGCAtgca]/; # This is the transliteration operator.
                    }
                    
                    
                    # Export the reference and query sequences
                    #-----------------------------------------
                    # Extract the query ID and sequence from the original FASTA.
                    my $query_sequence = $record[1]; # The first line in the FASTA is the sequence.
                    chomp $query_sequence;
                    print $paired_filehandle "@\n>$qseqid\n$query_sequence\n$reference_ID\n$reference_sequence\n"; # Output the query and the reference. Add the '>' back before the query ID and a '@' separator at the end.
                    }
		if ($batches > 0) {
			if ($batchround == 1) { # we have a batch coming through
				   @reference = split ("\>", $reference);
				my $entrycounter = 0;
				#my $entry = '';
				shift @reference;
				foreach $entry	(@reference) {
					my @entry = ();	
					@entry = split ("\n", $entry);
					$reference_ID = ">$entry[0]";
					$reference_sequence = join ('', @entry[1..$#entry]);
					unless ($blast_db ne '') {
						# efetch retrieved batches of sequences must then be excised from the whole entry
						#print "$reference_sequence, $fstarts[$entrycounter], $fstops[$entrycounter]\n";
						my @police_seq = (); @police_seq = split ('', $reference_sequence); my $police_size = ();
						$police_size = scalar @police_seq; my $police_stop = (); 
						if ($fstops[$entrycounter] > $police_size) { # the 3' end is out of range of the reference, just use as much sequence as possible
						#print "$reference_ID\n$reference_sequence, $fstarts[$entrycounter], $fstops[$entrycounter]\n";
						#exit;
						$police_stop = $police_size - 1;
						} else { 
						$police_stop = $fstops[$entrycounter];
						}
						$reference_sequence = excise_range ($reference_sequence, $fstarts[$entrycounter], $police_stop);
						#exit;
					}
					my $rev = '';
					if ($blast_db ne '') {
					  $rev = $reversed[$entrycounter];
					} else {
					  $rev = $freversed[$entrycounter];
					}
					if ($rev == 1) {
		                       		$reference_sequence = reverse($reference_sequence);
                			        $reference_sequence =~ tr/[ACGTacgt]/[TGCAtgca]/; # This is the transliteration operator.
					}
				if ($blast_db ne '') {
				$qseq = $query_seqs[$entrycounter];	
				$qseqid = $qseqid[$entrycounter];
				} else {
				# pull from the efetch prepared records which have large entries taken out at this stage
				$qseq = $fquery_seqs[$entrycounter];
				$qseqid = $fqseqids[$entrycounter];
				}
				chomp $qseq;
				print $paired_filehandle "@\n>$qseqid\n$qseq\n$reference_ID\n$reference_sequence\n"; # Output the query and the reference. Add the '>' back before the query ID and a '@' separator at the end.
				$entrycounter++;
				}
			# reset variables for next batch round
			$batchround = 0; @qseqid = (); @query_seqs = (); @reversed = (); @fqseqids = (); @fquery_seqs = (); @fstarts = (); @fstops = (); @freversed = (); @fsseqids = ();		
			}
		}
                } # The matching BLAST entry.
        } # while loop through the FASTA file.
	# if in batches, then it is likely there is a last batch to be processes with fewer entries than the defined batch size

               if ($batches > 0) {
			# need to check if there is a lagging batch to finish.
				my $qpile = '';
			if ($blast_db ne '') {
					$qpile = scalar @qseqid;
			} else {
					$qpile = scalar @fqseqids
			}

			if ($qpile > 0) {
			close BATCH; $batch = 0;
                        # pack up and ship out
		   if ($blast_db ne '') {	
                        $reference = `blastdbcmd -db $blast_db -entry_batch batchfile.$batchcount`;

                        if ($reference eq '') {
                        print "\n\nERROR: blastdbcmd cannot retrieve batchfile.$batchcount\n"; exit;
                        }
                       unlink ("batchfile.$batchcount");
		} else {

		
                        $fsseqids = join ("\,", @fsseqids);
                                        while (1) { # Run this loop until "last" is called.
#                                                $reference = `efetch -db Nucleotide -id $fsseqids -format fasta 2>> $error_log_filename`; # The efetch command.
                                                #print "\n\nEFETCH COMMAND: efetch -db nuccore -id $fsseqids -format fasta 2>> $error_log_filename";
                                                $reference = `efetch -db nuccore -id $fsseqids -format fasta 2>> $error_log_filename`; # The efetch command.

                                                if ($reference eq '') { # Failed efetches return a newline. If the return is completely empty, then efetch couldn't start for some reason - suggesting a different error.
                                                if ($first_efetch_start == 1) {
                                                print "\t\tERROR: efetch could not start. See $error_log_filename.\n";
                                                exit;
                                                } else {
                                                print "\t\tERROR: efetch could not start for sequence $qseqid. See $error_log_filename.\n";
                                                exit;
                                                }
                                        }

                                my @reference = split("\n", $reference);
                                if ( ($reference[1])  && ($reference =~ /\S/) ) { # If the reference was successfully split and contains at least one non-whitespace character, exit the loop.
                                $first_efetch_start = 0;
                                last;
                                }
                                }

		}			
                                @reference = split ("\>", $reference);
				shift @reference;
                                my $entrycounter = 0;
                                foreach $entry  (@reference) {
                                        my @entry = ();
                                        @entry = split ("\n", $entry);
                                        $reference_ID = ">$entry[0]";
                                        $reference_sequence = join ('', @entry[1..$#entry]);
					unless ($blast_db ne '') {
						# efetch retrieved batches of sequences must then be excised from the whole entry
						my @police_seq = (); @police_seq = split ('', $reference_sequence); my $police_size = ();
                                                $police_size = scalar @police_seq; my $police_stop = ();
                                                if ($fstops[$entrycounter] > $police_size) { # the 3' end is out of range of the reference, just use as much sequence as possible
                                                #print "$reference_ID\n$reference_sequence, $fstarts[$entrycounter], $fstops[$entrycounter]\n";
                                                #exit;
                                                $police_stop = $police_size - 1;
                                                } else {
                                                $police_stop = $fstops[$entrycounter];
                                                }
                                                $reference_sequence = excise_range ($reference_sequence, $fstarts[$entrycounter], $police_stop);

					}

					my $rev = '';
					if ($blast_db ne '') {
					  $rev = $reversed[$entrycounter];
					} else {
					  $rev = $freversed[$entrycounter];
					}
                                        if ($rev == 1) {
                                                $reference_sequence = reverse($reference_sequence);
                                                $reference_sequence =~ tr/[ACGTacgt]/[TGCAtgca]/; # This is the transliteration operator.
                                        }
				if ($blast_db ne '') {
                                $qseq = $query_seqs[$entrycounter];
                                $qseqid = $qseqid[$entrycounter];
				} else {
				# pull from the efetch prepared records which have large entries taken out at this stage
				$qseq = $fquery_seqs[$entrycounter];
				$qseqid = $fqseqids[$entrycounter];

				}
				chomp $qseq;
                                print $paired_filehandle "@\n>$qseqid\n$qseq\n$reference_ID\n$reference_sequence\n"; # Output the query and the reference. Add the '>' back before the query ID and a '@' separator at the end.
                                $entrycounter++;
                                }
		}
                

		# finally, if efetch option is invoked then we need to clear through large entries one by painful one
       
		unless ($blast_db ne '') {
			my $big = ''; my $bigcounter = 0; my $bigstart = ''; my $bigstop = ''; my $bigfqseqid = '';
			my $bigfquery_seq = ''; my $bigfreversed = '';
			my $checkcounter = 0;
			#foreach $big (@bigfsseqids) {print "bigfsseqid: $big bigfstart: $bigfstarts[$checkcounter] bigfstop: $bigfstops[$checkcounter] bigfqseqid: $bigfqseqids[$checkcounter]\n"; $checkcounter++;} exit; 
			if ($bigcheck == 1) {
				my $bigs = scalar @bigfsseqids;
				print "\tProcessing $bigs large entries....\n";
				foreach $big (@bigfsseqids) {
					$bigstart = $bigfstarts[$bigcounter];
					$bigstop = $bigfstops[$bigcounter];
					$bigfreversed = $bigfreversed[$bigcounter];
                    			while (1) { # Run this loop until "last" is called.
#                        			$reference = `efetch -db Nucleotide -id $big -seq_start $bigstart -seq_stop $bigstop -format fasta 2>> $error_log_filename`; # The efetch command.
                                                #print "\n\nEFETCH COMMAND: efetch -db nuccore -id $big -seq_start $bigstart -seq_stop $bigstop -format fasta 2>> $error_log_filename";
                                                $reference = `efetch -db nuccore -id $big -seq_start $bigstart -seq_stop $bigstop -format fasta 2>> $error_log_filename`; # The efetch command.

                        			if ($reference eq '') { # Failed efetches return a newline. If the return is completely empty, then efetch couldn't start for some reason - suggesting a different error.
                            			if ($first_efetch_start == 1) {
                                		print "\t\tERROR: efetch could not start. See $error_log_filename.\n";
                                		exit;
                            			} else {
                                		print "\t\tERROR: efetch could not start for sequence $qseqid. See $error_log_filename.\n";
                                		exit;
                            			}
                        		}

                        	my @reference = split("\n", $reference);
                        	if ( ($reference[1])  && ($reference =~ /\S/) ) { # If the reference was successfully split and contains at least one non-whitespace character, exit the loop.
                            	$first_efetch_start = 0;
                            	last;
                        	}
                    		}				 	
                    # Separate the header and put all sequence lines together.
                    @reference = split("\n", $reference);
                    $reference_ID = $reference[0];
                    $reference_sequence = join('', @reference[1..$#reference]);
    
                    # Is the reference section 3'-5'? If so, reverse complement it.
                    if ($bigfreversed == 1) {
                        $reference_sequence = reverse($reference_sequence);
                        $reference_sequence =~ tr/[ACGTacgt]/[TGCAtgca]/; # This is the transliteration operator.
                    }
                    
                    
                    # Export the reference and query sequences
                    #-----------------------------------------
                    # Extract the query ID and sequence from the original FASTA.
                    $bigfquery_seq = $bigfquery_seqs[$bigcounter]; # The first line in the FASTA is the sequence.
                    chomp $bigfquery_seq;
		    $bigfqseqid = $bigfqseqids[$bigcounter];	
		
                    print $paired_filehandle "@\n>$bigfqseqid\n$bigfquery_seq\n$reference_ID\n$reference_sequence\n"; # Output the query and the reference. Add the '>' back before the query ID and a '@' separator at the end.
 					


			$bigcounter++;	
				}
			}

		}
		}
 
        close $fasta_filehandle;
        $/ = "\n"; # Set the record separator back to the default newline.
        
        if ($count_reads_discarded > 0 ) {
            print "\t\t$count_reads_discarded read(s) excluded.\n";
        }
        
    } # Retrieving reference files.
    close $paired_filehandle; # We've finished writing to this, so close it. Then, when we open it below, we will start from the beginning of the file.
    
    
    ############### Step 3: Align queries and references #############  
    #my $paired_filename = "$basename.paired.txt"; # We already made this file in the last step.
    my $aln_filename = $output_directory . '/' .$basename . '.aln.txt';
    printf "\tAligning queries and references...\n";
     
    $/ = "@\n"; # Set the record separator to "@\n", which separates paired query-and-reference sequences in the input FASTA.
    
    open($paired_filehandle, $paired_filename) or die "\t\tCould not open paired file $paired_filename: $!\n";
    open(my $aln_filehandle, '>', $aln_filename) or die "\t\tCould not open $aln_filename for writing: $!\n";
    
    readline $paired_filehandle; # Skip the first record (just a "@").
    
    while (1) { # Run this loop until "last" is called.
                
        my $record = readline($paired_filehandle);
        if (! defined $record) { last }; # If the record is just a newline, exit the loop. You've processed the whole file.
    
        my @record = split("\n", $record); # Split the record into its four lines (the last element is just a newline; absent from the last record).
    
        my ($reference_aligned, $query_aligned) = needleman_wunsch($record[3], $record[1]); # Run the Needleman-Wunsch alignment algorithm. $record[3] is the reference and $record[1] is the query.
    
        print $aln_filehandle "$record[0]\t$query_aligned\t$reference_aligned\n"; # Print to the output file. $record[0] is the query name and base counts.
    }
    
    close $paired_filehandle;
    close $aln_filehandle;


 

    ############### Step 4: Summarise mismatches #############  
    my $alignment_count = 0; # Number of aligned sequences.
    my $C_count_at_position0 = 0;
    my $G_count_at_position0 = 0;
    
    open($aln_filehandle, $aln_filename) or die "\t\tCould not open alignment file $aln_filename: $!\n";
	$/ = "\n"; # reset record separator to newline bug fix
    my $mismatches_filename = $output_directory . '/' .$basename . '.mismatches.txt'; # Set up the output file.
    open (my $mismatches_filehandle, '>', $mismatches_filename) or die "Could not open $mismatches_filename for writing: $!\n";
    
    # Make empty arrays to store original and mismatch values for the 5'-3' direction and 3'-5' direction.
    my @As_5 = (0)x25;
    my @AtoTs_5 = (0)x25;
    my @AtoCs_5 = (0)x25;
    my @AtoGs_5 = (0)x25;
    my @Ts_5 = (0)x25;
    my @TtoAs_5 = (0)x25;
    my @TtoCs_5 = (0)x25;
    my @TtoGs_5 = (0)x25;
    my @Cs_5 = (0)x25;
    my @CtoAs_5 = (0)x25;
    my @CtoTs_5 = (0)x25;
    my @CtoGs_5 = (0)x25;
    my @Gs_5 = (0)x25;
    my @GtoAs_5 = (0)x25;
    my @GtoTs_5 = (0)x25;
    my @GtoCs_5 = (0)x25;
    
    my @As_3 = (0)x25;
    my @AtoTs_3 = (0)x25;
    my @AtoCs_3 = (0)x25;
    my @AtoGs_3 = (0)x25;
    my @Ts_3 = (0)x25;
    my @TtoAs_3 = (0)x25;
    my @TtoCs_3 = (0)x25;
    my @TtoGs_3 = (0)x25;
    my @Cs_3 = (0)x25;
    my @CtoAs_3 = (0)x25;
    my @CtoTs_3 = (0)x25;
    my @CtoGs_3 = (0)x25;
    my @Gs_3 = (0)x25;
    my @GtoAs_3 = (0)x25;
    my @GtoTs_3 = (0)x25;
    my @GtoCs_3 = (0)x25;
 
    # Search the alignment for any C->T and G->A mismatches on the 5' and 3' ends respectively.
    my $empty_aln_file = 1; # Flag. Default to empty.
    while (1) { # Run this loop until "last" is called.
        my $line = readline($aln_filehandle);

        if (!defined $line) { last; } # If there is no next line, exit.

        chomp $line;
        my @line = split("\t", $line); # Split the line into its four tab-delimited fields.
        
        if ($line[1] eq '') { # If there is no alignment:
            last; # The file is empty, so exit the while loop.
        } else {
            $empty_aln_file = 0; # This file isn't empty.
            $alignment_count ++; # Increment by 1.
            # Check for mismatches on the 5'-3' end:
            my @query_seq = split('', $line[1]); # Query sequence. Split into characters. Remember that the sequence runs from 5'-3'.
            my @subject_seq = split('', $line[2]); # Reference sequence.
            
            if ($subject_seq[0] eq 'C') { # Note if there's a C at position 0.
                $C_count_at_position0 = $C_count_at_position0 + 1;
            }
            
            foreach my $i (0..24) { # Run this loop 25 times: look at the first 25 bases in the alignment. Read the base in the reference and note any mismatches in the query.
                if ($subject_seq[$i] eq 'A') {
                    $As_5[$i]++; # If this base has a A in the subject, increase the value at this location in @As_5 by 1. Remember that all sequences can contribute to this total.
                    if ($query_seq[$i] eq 'T') { $AtoTs_5[$i]++; } # If this base also has a T in the query, increase the value at this location in @AtoTs_5 by 1. Remember that all sequences can contribute to this total.
                    if ($query_seq[$i] eq 'C') { $AtoCs_5[$i]++; } # If this base also has a C in the query, increase the value at this location in @AtoCs_5 by 1.
                    if ($query_seq[$i] eq 'G') { $AtoGs_5[$i]++; } # And so on...
                } elsif ($subject_seq[$i] eq 'T') {
                    $Ts_5[$i]++; # Increase Ts at this location by 1.
                    if ($query_seq[$i] eq 'A') { $TtoAs_5[$i]++; } # Note if there's a T to A mismatch. And so on.
                    if ($query_seq[$i] eq 'C') { $TtoCs_5[$i]++; }
                    if ($query_seq[$i] eq 'G') { $TtoGs_5[$i]++; }
                } elsif ($subject_seq[$i] eq 'C') {
                    $Cs_5[$i]++;
                    if ($query_seq[$i] eq 'A') { $CtoAs_5[$i]++; }
                    if ($query_seq[$i] eq 'T') { $CtoTs_5[$i]++; }
                    if ($query_seq[$i] eq 'G') { $CtoGs_5[$i]++; }
                } elsif ($subject_seq[$i] eq 'G') {
                    $Gs_5[$i]++;
                    if ($query_seq[$i] eq 'A') { $GtoAs_5[$i]++; }
                    if ($query_seq[$i] eq 'T') { $GtoCs_5[$i]++; }
                    if ($query_seq[$i] eq 'C') { $GtoCs_5[$i]++; }
                }
            } # 0-24 bases.
    
                
            # Check for mismatches on the 3'-5' end:
            my @query_seq_rev = split('', reverse($line[1]) ); # This time, reverse the sequences, so the first 25 positions are from the 3' end.
            my @subject_seq_rev = split('', reverse($line[2]) );
            if ($single_s == 0) {
            if ($subject_seq[0] eq 'G') { # Note if there's a G at position 0.
                $G_count_at_position0 = $G_count_at_position0 + 1;
            }
	    } else {
 	    
		if ($subject_seq[0] eq 'C') { # for single stranded preps check for C at position 0
		$G_count_at_position0 = $G_count_at_position0 + 1; #using the same variable as for G for double stranded libraries to avoid confusion with 5'C's
		}
	    }       
            foreach my $i (0..24) {
                if ($subject_seq_rev[$i] eq 'A') {
                    $As_3[$i]++;
                    if ($query_seq_rev[$i] eq 'T') { $AtoTs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'C') { $AtoCs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'G') { $AtoGs_3[$i]++; }
                } elsif ($subject_seq_rev[$i] eq 'T') {
                    $Ts_3[$i]++;
                    if ($query_seq_rev[$i] eq 'A') { $TtoAs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'C') { $TtoCs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'G') { $TtoGs_3[$i]++; }
                } elsif ($subject_seq_rev[$i] eq 'C') {
                    $Cs_3[$i]++; 
                    if ($query_seq_rev[$i] eq 'A') { $CtoAs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'T') { $CtoTs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'G') { $CtoGs_3[$i]++; }
                } elsif ($subject_seq_rev[$i] eq 'G') {
                    $Gs_3[$i]++; 
                    if ($query_seq_rev[$i] eq 'A') { $GtoAs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'T') { $GtoCs_3[$i]++; }
                    if ($query_seq_rev[$i] eq 'C') { $GtoCs_3[$i]++; }
                }
            } # 0-24 bases.
    
        } # If there was at least one alignment, search for mismatches.
    }

    if ($single_s == 0) {
    print $mismatches_filehandle "Number of aligned sequences: $alignment_count
Count of Cs at 5'-3' position 0: $C_count_at_position0
Count of Gs at 3'-5' position 0: $G_count_at_position0
Position\tP_{A->T}\tP_{A->C}\tP_{A->G}\tP_{T->A}\tP_{T->C}\tP_{T->G}\tP_{C->A}\tP_{C->T}\tP_{C->G}\tP_{G->A}\tP_{G->T}\tP_{G->C}\n"; # Only one column for each mismatch, because 5'-3' and 3'-5' will instead be differentiated by position. 3'-5' will have negative indices, like in mapDamage.;
} else {
    print $mismatches_filehandle "Number of aligned sequences: $alignment_count
Count of Cs at 5'-3' position 0: $C_count_at_position0
Count of Cs at 3'-5' position 0: $G_count_at_position0
Position\tP_{A->T}\tP_{A->C}\tP_{A->G}\tP_{T->A}\tP_{T->C}\tP_{T->G}\tP_{C->A}\tP_{C->T}\tP_{C->G}\tP_{G->A}\tP_{G->T}\tP_{G->C}\n"; # Only one column for each mismatch, because 5'-3' and 3'-5' will

}    
    if ($empty_aln_file == 1) {
        print "\t\t$aln_filename empty. Outputting NAs.\n";
        print $mismatches_filehandle "0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
1	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
2	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
3	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
4	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
5	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
6	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
7	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
8	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
9	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
10	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
11	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
12	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
13	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
14	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
15	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
16	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
17	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
18	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
19	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
20	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
21	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
22	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
23	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
24	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-24	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-23	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-22	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-21	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-20	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-19	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-18	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-17	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-16	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-15	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-14	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-13	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-12	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-11	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-10	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-9	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-8	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-7	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-6	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-5	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-4	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-3	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-2	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
-1	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
";
    } else {
        # Calculate and output mismatch scores at each base position. These variables will be overwritten for each position.
        # For mismatch XtoY, score = (number of XtoYs)/(number of Xs). Proportion of X bases in the subject that became Y in the query.
        my $mismatch_score_AtoT = 0;
        my $mismatch_score_AtoC = 0;
        my $mismatch_score_AtoG = 0;
        my $mismatch_score_CtoA = 0;
        my $mismatch_score_CtoT = 0;
        my $mismatch_score_CtoG = 0;
        my $mismatch_score_TtoA = 0;
        my $mismatch_score_TtoC = 0;
        my $mismatch_score_TtoG = 0;
        my $mismatch_score_GtoA = 0;
        my $mismatch_score_GtoT = 0;
        my $mismatch_score_GtoC = 0;
            
        # 5'-3' end:
        foreach my $i (0..24) {
            if ($As_5[$i] != 0) { # If this position had any As in the subject, calculate its mismatch scores.
                $mismatch_score_AtoT = $AtoTs_5[$i] / $As_5[$i];
                $mismatch_score_AtoC = $AtoCs_5[$i] / $As_5[$i];
                $mismatch_score_AtoG = $AtoGs_5[$i] / $As_5[$i];
            } if ($Ts_5[$i] != 0) {
                $mismatch_score_TtoA = $TtoAs_5[$i] / $Ts_5[$i];
                $mismatch_score_TtoC = $TtoCs_5[$i] / $Ts_5[$i];
                $mismatch_score_TtoG = $TtoGs_5[$i] / $Ts_5[$i];
            } if ($Cs_5[$i] != 0) {
                $mismatch_score_CtoA = $CtoAs_5[$i] / $Cs_5[$i];
                $mismatch_score_CtoT = $CtoTs_5[$i] / $Cs_5[$i];
                $mismatch_score_CtoG = $CtoGs_5[$i] / $Cs_5[$i];
            } if ($Gs_5[$i] != 0) {
                $mismatch_score_GtoA = $GtoAs_5[$i] / $Gs_5[$i];
                $mismatch_score_GtoT = $GtoTs_5[$i] / $Gs_5[$i];
                $mismatch_score_GtoC = $GtoCs_5[$i] / $Gs_5[$i];
            }
            print $mismatches_filehandle "$i\t$mismatch_score_AtoT\t$mismatch_score_AtoC\t$mismatch_score_AtoG\t$mismatch_score_TtoA\t$mismatch_score_TtoC\t$mismatch_score_TtoG\t$mismatch_score_CtoA\t$mismatch_score_CtoT\t$mismatch_score_CtoG\t$mismatch_score_GtoA\t$mismatch_score_GtoT\t$mismatch_score_GtoC\n"; # Print the position followed by all scores. For 5'-3', the positions are simply 0-24.
        }
        
        
        # 3'-5' end (reset scores first):
        $mismatch_score_AtoT = 0;
        $mismatch_score_AtoC = 0;
        $mismatch_score_AtoG = 0;
        $mismatch_score_CtoA = 0;
        $mismatch_score_CtoT = 0;
        $mismatch_score_CtoG = 0;
        $mismatch_score_TtoA = 0;
        $mismatch_score_TtoC = 0;
        $mismatch_score_TtoG = 0;
        $mismatch_score_GtoA = 0;
        $mismatch_score_GtoT = 0;
        $mismatch_score_GtoC = 0;
    
        foreach my $i (reverse(0..24)) {
            if ($As_3[$i] != 0) { # If this position had any As, calculate its mismatch scores.
                $mismatch_score_AtoT = $AtoTs_3[$i] / $As_3[$i];
                $mismatch_score_AtoC = $AtoCs_3[$i] / $As_3[$i];
                $mismatch_score_AtoG = $AtoGs_3[$i] / $As_3[$i];
            } if ($Ts_3[$i] != 0) {
                $mismatch_score_TtoA = $TtoAs_3[$i] / $Ts_3[$i];
                $mismatch_score_TtoC = $TtoCs_3[$i] / $Ts_3[$i];
                $mismatch_score_TtoG = $TtoGs_3[$i] / $Ts_3[$i];
            } if ($Cs_3[$i] != 0) {
                $mismatch_score_CtoA = $CtoAs_3[$i] / $Cs_3[$i];
                $mismatch_score_CtoT = $CtoTs_3[$i] / $Cs_3[$i];
                $mismatch_score_CtoG = $CtoGs_3[$i] / $Cs_3[$i];
            } if ($Gs_3[$i] != 0) {
                $mismatch_score_GtoA = $GtoAs_3[$i] / $Gs_3[$i];
                $mismatch_score_GtoT = $GtoTs_3[$i] / $Gs_3[$i];
                $mismatch_score_GtoC = $GtoCs_3[$i] / $Gs_3[$i];
            }
            my $position = 0 - $i; # For 3'-5', make the position negative.
            print $mismatches_filehandle "$position\t$mismatch_score_AtoT\t$mismatch_score_AtoC\t$mismatch_score_AtoG\t$mismatch_score_TtoA\t$mismatch_score_TtoC\t$mismatch_score_TtoG\t$mismatch_score_CtoA\t$mismatch_score_CtoT\t$mismatch_score_CtoG\t$mismatch_score_GtoA\t$mismatch_score_GtoT\t$mismatch_score_GtoC\n"; # Print the position followed by all scores.
        }
    } # Print non-NA outputs.
    
    close $aln_filehandle;
    close $mismatches_filehandle;
    
    
    ############### Step 5: Plot mismatches ############# 
    print "\tPlotting mismatches...\n";
    my $r_script = '';
    if ($single_s == 0) { 
        #$r_script= "Plot_mismatches_ds.R";
	`./Plot_mismatches_ds.R $mismatches_filename`;
    } else {
	#$r_script= "Plot_mismatches_ss.R";
	`./Plot_mismatches_ss.R $mismatches_filename`;
    }
    #`Rscript $r_script $mismatches_filename`;    
    unlink "Rplots.pdf";
    print "\t################ MetaDamage analysis complete...Go to \\$output_directory for results ################\n"; 

$counter++;    
} # For each input FASTA.



#####################
#### SUBROUTINES ####
#####################
sub needleman_wunsch {
    # Adapted from "Needleman-Wunsch Algorithm", Copyright (C) 2013-2020 S. Evan Staton, https://github.com/sestaton/sesbio/blob/master/phylogenetics/needleman-wunsch.pl

    # usage statement
    #die "usage: $0 <sequence 1> <sequence 2>\n" unless @ARGV == 2;
    
    # get sequences from command line
    #my ($seq1, $seq2) = @ARGV;
    my ($seq1, $seq2) = @_;

    # scoring scheme
    my $MATCH    =  1; # +1 for letters that match
    my $MISMATCH = -1; # -1 for letters that mismatch
    my $GAP      = -1; # -1 for any gap
    
    # initialization
    my @matrix;
    $matrix[0][0]{score}   = 0;
    $matrix[0][0]{pointer} = "none";
    for (my $j = 1; $j <= length($seq1); $j++) {
        $matrix[0][$j]{score}   = $GAP * $j;
        $matrix[0][$j]{pointer} = "left";
    }
    for (my $i = 1; $i <= length($seq2); $i++) {
        $matrix[$i][0]{score}   = $GAP * $i;
        $matrix[$i][0]{pointer} = "up";
    }
    
    # fill
    for (my $i = 1; $i <= length($seq2); $i++) {
        for (my $j = 1; $j <= length($seq1); $j++) {
            my ($diagonal_score, $left_score, $up_score);
    
            # calculate match score
            my $letter1 = substr($seq1, $j-1, 1);
            my $letter2 = substr($seq2, $i-1, 1);                            
            if ($letter1 eq $letter2) {
                $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
            }
            else {
                $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
            }
    
            # calculate gap scores
            $up_score   = $matrix[$i-1][$j]{score} + $GAP;
            $left_score = $matrix[$i][$j-1]{score} + $GAP;
    
            # choose best score
            if ($diagonal_score >= $up_score) {
                if ($diagonal_score >= $left_score) {
                    $matrix[$i][$j]{score}   = $diagonal_score;
                    $matrix[$i][$j]{pointer} = "diagonal";
                }
            else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "left";
                }
            } 
            else {
                if ($up_score >= $left_score) {
                    $matrix[$i][$j]{score}   = $up_score;
                    $matrix[$i][$j]{pointer} = "up";
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "left";
                }
            }
        }
    }
    
    # trace-back
    
    my $align1 = "";
    my $align2 = "";
    
    # start at last cell of matrix
    my $j = length($seq1);
    my $i = length($seq2);
    
    while (1) {
        last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix
    
        if ($matrix[$i][$j]{pointer} eq "diagonal") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "left") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= "-";
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "up") {
            $align1 .= "-";
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
        }    
    }
    
    $align1 = reverse $align1;
    $align2 = reverse $align2;
    return ($align1, $align2);
}

sub excise_range {
		my $seq = $_[0];
		my $start = $_[1];
		my $stop = $_[2];
		my @seq = (); my $base = (); my $basecounter = 0; 
		my @excseq = (); my $length = $stop-$start+1;
		my $add = $start-1; my $excseq = ();
		@seq = split ('', $seq);
		while ($basecounter < $length) {
			push (@excseq, $seq[$add]);
			$add++; $basecounter++;
		}
		$excseq = join ('', @excseq);
	return $excseq;
}
