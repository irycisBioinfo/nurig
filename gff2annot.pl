#!/usr/bin/perl

use strict;
use Getopt::Long;

my @gff;
my @c;
my $l;
my $label;
my $gff;
my $outFile;
my @Options;
my $hypo;
my @fields;
my @c2;
my $term;
my $f;
my $type;


@Options = (
		
		{OPT=>"gff=s",	VAR=>\$gff,	DESC=>"GFF File"},
		{OPT=>"out=s",	VAR=>\$outFile,	DEFAULT => "GFF_File.tab", DESC=>"Output annot file"},
		{OPT=>"hypo!",	VAR=>\$hypo, DEFAULT => 1, DESC=>"Remove hypothetical proteins? (default TRUE)"},
		{OPT=>"fields=s{,}", VAR=>\@fields, DESC=>"Space-separate list of terms to extract from GFF attribute field (in order of preference)"},
		{OPT=>"type=s", VAR=>\$type, DEFAULT => "CDS", DESC => "Type of feature to parse (e.g CDS, gene, exon...) default: CDS)"},
		);


#Check options and set variables
(@ARGV < 1) && (usage());
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}

open(GFF,$gff);
@gff = <GFF>;
close GFF;


open(OUT,">$outFile");

foreach $l (@gff)
{
	if($l =~ /\t$type\t/)
	{
		$label = "";
		@c = split("\t",$l);
		
		
		
		@c2 = split(/;/,$c[-1]);
		LOOP:{
			
			foreach $f (@fields)
			{
				foreach $term (@c2)
				{
					if($term =~ /$f=/)
					{
						$label = $term;
						$label =~ s/$f=//;
						if($hypo & $label =~ /hypothetical protein/i)
						{
							$label = "";
						}
						
						last LOOP;
					}
				}
			}
		}
		print OUT "$c[0]\t$c[3]\t$c[4]\t$c[6]\t$label\n";
	}
}

sub usage(){
	print "Usage:\n\n";
	print "Gff2annot\n";
	print "writen by: Val F. Lanza (valfernandez.vf\@gmail.com)\n\n";
	
	print "\nParameters:\n\n";
	foreach (@Options) {
		
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
	
}
