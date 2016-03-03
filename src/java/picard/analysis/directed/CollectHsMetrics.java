/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Collects a set of HS metrics from a sam or bam file.  See HsMetricsCollector and CollectTargetedMetrics for more details.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = CollectHsMetrics.USAGE_SUMMARY + CollectHsMetrics.USAGE_DETAILS,
        usageShort = CollectHsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectHsMetrics extends CollectTargetedMetrics<HsMetrics, HsMetricCollector> {
        static final String USAGE_SUMMARY = "Collects hybrid-selection (HS) specific metrics for a SAM or BAM file.  ";
        static final String USAGE_DETAILS = "<p>These metrics enable users to determine the efficacy and quality of hybrid selection (HS) experiments.</p>  " +
                "" +
                "<p>Hybrid selection enables users to target and sequence specific regions of a genome.  It is commonly used to characterize exon sequences" +
                " from genomic DNA or filter out contaminating bacterial DNA sequences from clinical samples.  For additional information " +
                "and the theory behind this technique, please see the following GATK " +
                "<a href=\"http://www.broadinstitute.org/gatk/guide/article?id=6331\">Dictionary</a> entry.</p>" +
                "" +
                "<p>This tool requires an aligned SAM or BAM file as well as bait and target interval files.  These interval files can be created from BED files using Picard's " +
                "<a href=\"http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList\">BedToInterval</a> tool.</p>" +
                "" +
                "If a reference sequence is provided, this program will calculate both AT_DROPOUT and GC_DROPOUT metrics.  These \"dropout\" " +
                "metrics are an attempt to measure the reduced representation of reads, which contain sequences that deviate from 50% G/C content.  " +
                "This reduction in the number of aligned reads is due to the increased numbers of errors associated with sequencing " +
                "regions of excessive or deficient numbers of G/C bases, ultimately leading to poor mapping efficiencies." +
                "" +
                "<p>The PER_TARGET_COVERAGE option can be invoked to output G/C content and mean sequence depth information for every target interval." +
                "</p>" +
                "" +
                "<p>Please note that coverage measurements are capped at ~32K to constrain memory usage."+

                "<h4>Usage Example:</h4>"+
                "<pre>" +
                "java -jar picard.jar CollectHsMetrics \\<br />" +
                "      I=input.bam \\<br />" +
                "      O=hs_metrics.txt \\<br />" +
                "      R=reference_sequence.fasta \\<br />" +
                "      BAIT_INTERVALS=bait.interval_list \\<br />" +
                "      TARGET_INTERVALS=target.interval_list" +
                "</pre> "   +
                "<p>Please see the CollectTargetedMetrics " +
                "<a href=\"http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics\">definitions</a> for a " +
                "complete description of the metrics produced by this tool.</p> "+
                "<hr />";

    @Option(shortName = "BI", doc = "An interval list file that contains the locations of the baits used.", minElements=1)
    public List<File> BAIT_INTERVALS;

    @Option(shortName = "N", doc = "Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional = true)
    public String BAIT_SET_NAME;

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(doc = "True if we are to clip overlapping reads, false otherwise.", optional=true, overridable = true)
    public boolean CLIP_OVERLAPPING_READS = true;

    @Override
    protected IntervalList getProbeIntervals() {
        for (final File file : BAIT_INTERVALS) IOUtil.assertFileIsReadable(file);
        return IntervalList.fromFiles(BAIT_INTERVALS);
    }

    @Override
    protected String getProbeSetName() {
        if (BAIT_SET_NAME != null) {
            return BAIT_SET_NAME;
        } else {
            final SortedSet<String> baitSetNames = new TreeSet<String>();
            for (final File file : BAIT_INTERVALS) {
                baitSetNames.add(CollectTargetedMetrics.renderProbeNameFromFile(file));
            }
            return StringUtil.join(".", baitSetNames);
        }
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    @Override
    protected HsMetricCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                              final List<SAMReadGroupRecord> samRgRecords,
                                              final ReferenceSequenceFile refFile,
                                              final File perTargetCoverage,
                                              final File perBaseCoverage,
                                              final IntervalList targetIntervals,
                                              final IntervalList probeIntervals,
                                              final String probeSetName,
                                              final int nearProbeDistance) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance,
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, CLIP_OVERLAPPING_READS, true, COVERAGE_CAP, SAMPLE_SIZE);
    }
}