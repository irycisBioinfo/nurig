#!/usr/bin/perl

use Cwd 'abs_path';
use Statistics::R;
use strict;
use Getopt::Long;

my $ref;
my $outFile;
my $table;
my $palette;
my $W;
my $H;
my $format;
my $PATH;
my $R;
my $cmdLista;
my @Options;

@Options = (
		
		{OPT=>"ref=s",	VAR=>\$ref,	DESC=>"Reference Genome"},
		{OPT=>"out=s",	VAR=>\$outFile,	DEFAULT => 'out', DESC=>"Network filename"},
		{OPT=>"table=s",	VAR=>\$table, DESC=>"Table with genome information"},
		{OPT=>"palette=s",	VAR=>\$palette,	DEFAULT => 'rainbow', DESC=>"Color palette for figure (rainbow, topo, terrain)"},
		{OPT=>"W=i",	VAR=>\$W,	DEFAULT => '1000', DESC=>"Image weight"},
		{OPT=>"H=i",	VAR=>\$H,	DEFAULT => '1000', DESC=>"Image height"},
		{OPT=>"format=s",	VAR=>\$format,	DEFAULT => 'SVG', DESC=>"Output format image (JPG,PNG,SVG,SVGZ)"}
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
	

unless($ref){
	print STDERR "You must specified a reference genome file\n";
	&usage();
}

unless($table){
	print STDERR "You must specified table genomes file\n";
	&usage();
}

$PATH = abs_path($0);
$PATH =~ s/\/NuRIG.pl//;

$R = Statistics::R->new(shared => 1);
$R->start();
$R->set('ref',$ref);
$cmdLista = "lista = as.data.frame(read.table('$table'))";
$R->run($cmdLista);
$R->set('palette',$palette);

$R->run_from_file("$PATH/lib/NuRIG.R");
$R-> stop();

system("java -jar $PATH/bin/cgview/cgview.jar -i out.xml -o $outFile.$format -f $format -W $W -H $H");



sub usage {
	print "Usage:\n\n";
	print "Nucmer Image Ring Generator (NuRIG) BETA-Version\n";
	print "writen by: Val F. Lanza (valfernandez.vf\@gmail.com)\n\n";
	print "NuRIG.pl --ref reference.fasta --table list.txt --out image.svg --palette rainbow --W 1000 --H 1000 --format SVG\n";

	print "\nParameters:\n\n";
	foreach (@Options) {
		
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
}
