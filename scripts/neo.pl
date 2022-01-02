#!/usr/bin/perl 

# This script reads a list of junctions from $ARGV[0] and parses 
# a SAM from STDIN to extract and print neojunctions in STAR chimeric 
# format: chr1 pos1 str1 chr2 pos2 str2 read_id
# if $ARGV[1] is 1 then the strand is additionally flipped
# P.S. is it better to flip mate1 instead of mate2?


$BAM_FREAD1 = 0x40;
$BAM_FREAD2 = 0x80;
$BAM_REVERSE = 0x10;
$BAM_SECONDARY = 0x100; 

open FILE, $ARGV[0];
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end) = split /\t/, $line;
    $intron{$chr}{$beg}{$end}++;
}
close FILE;

$mate = $ARGV[1];


@STRAND = ("+","-");
while(<STDIN>) {
    ($id, $flag, $ref, $pos, $qual, $cigar) = split /\t/;
    next if($flag & $BAM_SECONDARY); # skip secondary alignments

    $rev = ($flag & $BAM_REVERSE) ? 1 : 0; # read reverse complemented yes no
    $str = $STRAND[($rev + $mate + 1) & 1]; # additionally reverse complement if mate0 

    while($cigar=~/(\d+)(\w)/g) {
        $increment = $1;
        $operation = $2;
        if($operation eq 'M') {  
            $pos += $increment;
        }
        if($operation eq 'D') {  
            $pos += $increment;
        }
        if($operation eq 'N') {
            $beg = $pos;
            $end = $pos + $increment - 1;
	    unless($intron{$ref}{$beg}{$end}) {
		($beg, $end) = reverse ($beg, $end) if($str eq "-");
            	print join("\t", $ref, $beg, $beg + 1, $str, $ref, $end, $end + 1, $str, $id, $mate), "\n";
	    }
            $pos += $increment; 
        }
    }
}
