#!/usr/bin/perl

use Cwd 'abs_path';
use Statistics::R;

$PATH = abs_path($0);
$PATH =~ s/\/NuRIGWeb.pl//;

$R = Statistics::R->new(shared => 1);
$R->send("setwd('$PATH')");
$R->send("shiny::runApp('lib/app.r', launch.browser = TRUE)");
$R-> stop();


