#!/usr/bin/env perl
use strict;
use Getopt::Long;
use FileHandle;
use Carp;
use Data::Dumper;
use Storable;
use Storable qw(dclone);
use Text::NSP::Measures::2D::Fisher::twotailed ();
use Text::NSP::Measures::2D::Fisher::left      ();
use Text::NSP::Measures::2D::Fisher::right     ();

my ( $infile, $outfile, $p_cutoff );
GetOptions(
    "i=s"        => \$infile,
    "o=s"        => \$outfile,
    "p_cutoff=f" => \$p_cutoff,
);

#QA
croak "FATAL: you must provide input file -$infile-"
  unless $infile && -e $infile;
$p_cutoff = 0.05                          unless defined($p_cutoff);
$outfile  = $infile . '_ReEval_NegStatus' unless defined($outfile);
croak "FATAL: p_cutoff $p_cutoff incorrect"
  unless ( 0 < $p_cutoff and $p_cutoff <= 1 );

#load input
my $r;
my $r_in = w_read($infile);
my ( $id_header, @samples ) = split( "\t", $r_in->[0] );
my $N_samples_header = scalar @samples;
for my $i ( 1 .. $#{$r_in} ) {
    my $line = $r_in->[$i];
    my ( $id, @cols ) = split( "\t", $line );
    my $N_cols = scalar @cols;
    croak
"FATAL: column number in line $line $N_cols not match to header $N_samples_header"
      unless $#cols == $#samples;
    for my $j ( 0 .. $#cols ) {
        my $sample          = $samples[$j];
        my $status_combined = $cols[$j];
        my ( $mut, $ref, $status ) = split( "\,", $status_combined );

        #QA
        croak
          "mut read count is not an integer in status_combined $status_combined"
          unless $mut =~ /^\d+$/;
        croak
          "ref read count is not an integer in status_combined $status_combined"
          unless $ref =~ /^\d+$/;
        $r->{$id}{'insert_counts'}{$sample}{'mut'} = $mut;
        $r->{$id}{'insert_counts'}{$sample}{'ref'} = $ref;
        if ( defined($status) and $status eq '+' ) {
            $r->{$id}{'insert_counts'}{$sample}{'status'} = 1;
        }
        else {
            if ( $mut >= 3 ) {
                if ( $mut / ( $mut + $ref ) > 0.2 ) {
                    warn
"WARN: mutation $id sample $sample possibly positive but not labeled as positive mut $mut ref $ref\n";
                }
            }
        }
        $r->{$id}{'insert_counts'}{$sample}{'status'} = 1
          if ( defined($status) and $status eq '+' );
    }
}
re_evaluate_negative_samples_by_FET( $r, { 'cutoff' => $p_cutoff } );

#output
my @out = ( $r_in->[0] );
for my $i ( 1 .. $#{$r_in} ) {
    my $line = $r_in->[$i];
    my ( $id, @cols ) = split( "\t", $line );
    my @t = ($id);
    for my $j ( 0 .. $#cols ) {
        my $sample = $samples[$j];
        my $mut =
          deep_defined( $r, $id, 'insert_counts', $sample, 'mut' )
          ? $r->{$id}{'insert_counts'}{$sample}{'mut'}
          : ( print Dumper $r->{$id}{'insert_counts'}
              and croak "FATAL:cannot find mut for id $id sample $sample" );
        my $ref =
          deep_defined( $r, $id, 'insert_counts', $sample, 'ref' )
          ? $r->{$id}{'insert_counts'}{$sample}{'ref'}
          : ( print Dumper $r->{$id}{'insert_counts'}
              and croak "FATAL:cannot find ref for id $id sample $sample" );
        my $status =
          deep_defined( $r, $id, 'insert_counts', $sample, 'status' )
          ? $r->{$id}{'insert_counts'}{$sample}{'status'}
          : ( print Dumper $r->{$id}{'insert_counts'}
              and croak "FATAL:cannot find status for id $id sample $sample" );
        $status =~ s/^-$/Unknown/;
        $status =~ s/^0$/\-/;
        $status =~ s/^1$/\+/;
        push( @t, join( ",", ( $mut, $ref, $status ) ) );
    }
    push( @out, join( "\t", @t ) );
}

#output
w_write( \@out, $outfile, 1 );

sub w_read {
    my ($infile) = @_;
    croak "cannot find infile $infile" unless ( $infile && -e $infile );
    my @out;
    my $fh_in = FileHandle->new($infile) or croak "cannot open $infile";
    while (<$fh_in>) {
        chomp;
        push( @out, $_ );
    }
    $fh_in->close;
    my $in_ref = \@out;
    return $in_ref;
}

sub w_write {
    my ( $out_ref, $outfile, $option ) = @_;
    $option = 0 unless $option;
    croak "not an array ref, cannot write" unless ref($out_ref) eq 'ARRAY';
    croak "you must specify an output file to write" unless $outfile;
    my $bk_file = $outfile . "_bk_" . RandomString();
    my $md5_old = '';
    croak "bk_file already exist, unlikely with random string $bk_file"
      if -e $bk_file;
    if ( -e $outfile ) {
        croak "output file already exist, cannot use option 0" if $option == 0;
        if ( $option == 3 || $option == 4 ) {
            $md5_old = md5_file($outfile);
            rename $outfile, $bk_file;
        }
    }
    my $fh_out;
    if ( $option == 0 || $option == 1 || $option == 3 || $option == 4 ) {
        $fh_out = FileHandle->new("> $outfile")
          or croak "cannot write to output file $outfile";
    }
    elsif ( $option == 2 ) {
        $fh_out = FileHandle->new(">> $outfile")
          or croak "cannot write to output file $outfile";
    }
    else {
        croak "unrecognized option $option";
    }
    for ( @{$out_ref} ) {
        chomp;
        print $fh_out "$_\n";
    }
    if ( -e $bk_file ) {
        my $md5_new = md5_file($outfile);
        if ( $md5_old eq $md5_new ) {
            if ( $option == 3 ) {
                print
"output opt #3 - keep older output, remove newer outfile with the same md5\n";
                unlink($outfile);
                rename $bk_file, $outfile;
            }
            elsif ( $option == 4 ) {
                print
"output opt #4 - keep newer output, remove older outfile with the same md5\n";
                unlink($bk_file);
            }
            else {
                croak "FATAL: Strange error, option is not 3 or 4, -$option-";
            }
        }
    }
    $fh_out->close;
}

sub FisherExact
{ #http://search.cpan.org/dist/Text-NSP/lib/Text/NSP/Measures/2D/Fisher/twotailed.pm
    croak "input are not 5 variables" unless ( scalar @_ == 5 );
    my $npp = $_[0] + $_[1] + $_[2] + $_[3];
    return 1 if $npp == 0;    #wl modified 2011_05_06
    my $n1p  = $_[0] + $_[1];
    my $np1  = $_[0] + $_[2];
    my $n11  = $_[0];
    my $type = $_[4];
    my $Pvalue;

    if ( $type eq 'left' ) {
        $Pvalue = Text::NSP::Measures::2D::Fisher::left::calculateStatistic(
            n11 => $n11,
            n1p => $n1p,
            np1 => $np1,
            npp => $npp
        );
    }
    elsif ( $type eq 'right' ) {
        $Pvalue = Text::NSP::Measures::2D::Fisher::right::calculateStatistic(
            n11 => $n11,
            n1p => $n1p,
            np1 => $np1,
            npp => $npp
        );
    }
    elsif ( $type eq 'twotailed' ) {
        $Pvalue =
          Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
            n11 => $n11,
            n1p => $n1p,
            np1 => $np1,
            npp => $npp
          );
    }
    else {
        croak "incorrect Fisher test type $type";
    }
    if (
        (
            my $errorCode =
            Text::NSP::Measures::2D::Fisher::twotailed::getErrorCode()
        )
      )
    {
        print STDERR $errorCode . " - "
          . Text::NSP::Measures::2D::Fisher::twotailed::getErrorMessage();
    }
    else {
        return $Pvalue;
    }

}

sub deep_defined {    #avoid autovivification
    my ( $ref, @keys ) = @_;
    unless (@keys) {
        warn "deep_defined: no keys";
        return;
    }
    foreach my $key (@keys) {
        return unless defined($key);
        if ( ref $ref eq 'HASH' ) {
            return unless defined( $ref->{$key} );
            $ref = $ref->{$key};
            next;
        }
        if ( ref $ref eq 'ARRAY' ) {
            return unless 0 <= $key && $key < @{$ref};
            return unless defined( $ref->[$key] );
            $ref = $ref->[$key];
            next;
        }
        return;
    }
    return 1;
}

sub re_evaluate_negative_samples_by_FET {
    my ( $r, $opts ) = @_;
    croak "Fatal: input must be hash ref" unless ref($r) eq 'HASH';
    my $cutoff = defined( $opts->{'cutoff'} ) ? $opts->{'cutoff'} : 0.05;
    for my $id ( keys %{$r} ) {
        my ( @samples_positive, @samples_nonpositive );
        for my $sample_t ( keys %{ $r->{$id}{'insert_counts'} } ) {
            if ( deep_defined( $r, $id, 'insert_counts', $sample_t, 'status' )
                and $r->{$id}{'insert_counts'}{$sample_t}{'status'} eq '1' )
            {
                push( @samples_positive, $sample_t );
            }
            else {
                push( @samples_nonpositive, $sample_t );
            }
        }
        my $n_positive = scalar @samples_positive;
        croak
          "FATAL: input mutation $id must be positive in at least one samples"
          unless $n_positive >= 1;
        for my $sample2 (@samples_nonpositive) {
            my $mut2 =
              deep_defined( $r, $id, 'insert_counts', $sample2, 'mut' )
              ? $r->{$id}{'insert_counts'}{$sample2}{'mut'}
              : croak "FATAL: failed t
o find mut for id $id in sample $sample2";
            my $ref2 =
              deep_defined( $r, $id, 'insert_counts', $sample2, 'ref' )
              ? $r->{$id}{'insert_counts'}{$sample2}{'ref'}
              : croak "FATAL: failed t
o find ref for id $id in sample $sample2";
          CHECK_POSITIVE: for my $sample1 (@samples_positive) {
                my $mut1 =
                  deep_defined( $r, $id, 'insert_counts', $sample1, 'mut' )
                  ? $r->{$id}{'insert_counts'}{$sample1}{'mut'}
                  : croak "FATAL:
failed to find mut for id $id in sample $sample1";
                my $ref1 =
                  deep_defined( $r, $id, 'insert_counts', $sample1, 'ref' )
                  ? $r->{$id}{'insert_counts'}{$sample1}{'ref'}
                  : croak "FATAL: failed to fi
nd ref for id $id in sample $sample1";
                my $p = FisherExact( $mut2, $ref2, $mut1, $ref1, 'twotailed' );
                my $cutoff2 =
                  defined( $opts->{'Bonferroni_correction'} )
                  ? $cutoff * $n_positive
                  : $cutoff;
                if ( defined( $opts->{'keep_raw_record'} ) ) {
                    $r->{$id}{'insert_counts'}{$sample2}{'raw_record'}{$sample1}
                      = {
                        'p'          => $p,
                        'n_positive' => $n_positive,
                        'cutoff2'    => $cutoff2,
                        'mut2-ref2-mut1-ref1' =>
                          join( "_", ( $mut2, $ref2, $mut1, $ref1 ) ),
                      };
                }
                unless ( $p < $cutoff2 ) {
                    $r->{$id}{'insert_counts'}{$sample2}{'status'} = '-';
                    last CHECK_POSITIVE;
                }
            }
            $r->{$id}{'insert_counts'}{$sample2}{'status'} = 0
              unless deep_defined( $r, $id, 'insert_counts', $sample2,
                'status' );
        }
    }
}

sub RandomString {
    my ($l) = @_;
    $l = 10 unless $l;
    return join( "", map { ( "a" .. "z" )[ rand 26 ] } ( 1 .. $l ) );
}
