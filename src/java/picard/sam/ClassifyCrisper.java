package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author kcibul@gmail.com
 */
@CommandLineProgramProperties(
        usage = CleanSam.USAGE,
        usageShort = CleanSam.USAGE,
        programGroup = SamOrBam.class
)
public class ClassifyCrisper extends CommandLineProgram {

    static final String USAGE = "Classify impact to CRISPER target region";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input BAM")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write reads post classification")
    public File OUTPUT;

    @Option(doc = "minimum base quality score to consider for a mismatch/mutation")
    public Integer MIN_BASE_QUAL = 20;

    @Option(doc = "contig of target")
    public String TARGET_CONTIG;

    @Option(doc = "start position of target (1-base, inclusive)" )
    public Integer TARGET_START;

    @Option(doc = "end position of target (1-base, inclusive)" )
    public Integer TARGET_END;

    @Option(doc = "only consider this read name, for debugging", optional = true)
    public String DEBUG_READ_NAME;

    final private Map<String,Integer> summaryCounts = new HashMap<>();

    public static void main(final String[] argv) {
        new CleanSam().instanceMainWithExit(argv);
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReaderFactory factory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);

        if (VALIDATION_STRINGENCY == ValidationStringency.STRICT) {
            factory.validationStringency(ValidationStringency.LENIENT);
        }
        final SamReader reader = factory.open(INPUT);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        final CloseableIterator<SAMRecord> it = reader.iterator();
        final ProgressLogger progress = new ProgressLogger(Log.getInstance(CleanSam.class));

        // If the read (or its mate) maps off the end of the alignment, clip it
        while (it.hasNext()) {
            final SAMRecord rec = it.next();

            // says 'clipped'... but doesn't look it...
//            String debugReadName = "HWI:1:X:1:2114:9900:10252";
            if (DEBUG_READ_NAME != null && !rec.getReadName().equals(DEBUG_READ_NAME)) {
                continue;
            }

            // If the read (or its mate) maps off the end of the alignment, clip it
            AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference(rec);

            // check the read's mapping quality
            if (rec.getReadUnmappedFlag() && 0 != rec.getMappingQuality()) {
                rec.setMappingQuality(0);
            }

            int targetLength = TARGET_END - TARGET_START + 1;

            String type;

            // check that the target region is completely within the read
            if (rec.getContig() != null && rec.getContig().equals(TARGET_CONTIG) && rec.getAlignmentStart() <= TARGET_START && rec.getAlignmentEnd() >= TARGET_END) {

                int start = rec.getReadPositionAtReferencePosition(TARGET_START);
                int end = rec.getReadPositionAtReferencePosition(TARGET_END);

                int currentPosition = rec.getAlignmentStart();
                int insertedBases = 0;
                int deletedBases = 0;

                for (CigarElement ce : rec.getCigar().getCigarElements()) {
                    if (ce.getOperator() == CigarOperator.HARD_CLIP | ce.getOperator() == CigarOperator.SOFT_CLIP) {
                        continue;
                    }

                    if (ce.getOperator() == CigarOperator.M) {
                        currentPosition += ce.getLength();
                    }

                    if (ce.getOperator() == CigarOperator.INSERTION) {
                        insertedBases += ce.getLength();
                        // don't increment position
                    }

                    if (ce.getOperator() == CigarOperator.DELETION) {
                        deletedBases += ce.getLength();
                        currentPosition += ce.getLength();
                    }

                }

                // we didn't make it through the targeted region
                if (currentPosition <= TARGET_END && insertedBases + deletedBases == 0) {
                    type = "no-overlap"; // would be strange...
                } else if (insertedBases == 0 && deletedBases == 0) {

                    if (start >= 0 && end >= start ) {

                        // start and end are 1-based, but we need them to be zero based here
                        String bases = rec.getReadString().substring(start-1, end-1+1);  // substring is end exclusive
                        byte[] quals = Arrays.copyOfRange(rec.getBaseQualities(), start-1, end-1+1);  // copyOfRange is end exclusive

                        int mutatedBases = 0;
                        int noisyBases = 0;

                        if (bases.length() < targetLength) {
                            System.out.println("Error!");
                        }
                        for (int i = 0; i < targetLength; i++) {
                            if ('=' != bases.charAt(i)) {
                                if (quals[i] > MIN_BASE_QUAL) {
                                    mutatedBases++;
                                } else {
                                    noisyBases++;
                                }
                            }
                        }

                        if (mutatedBases > 0) {
                            type = "mutation";
                        } else if (noisyBases > 0) {
                            type = "WT-noise";
                        } else {
                            type = "WT";
                        }
                    } else {
                        type = "clipped";
                    }
                } else if (insertedBases - deletedBases % 3 == 0) {
                    type = "in-frame";
                } else {
                    type = "frame-shift";
                }


            } else {
                // if we are within 250bp, call it 'near-target'
                if (rec.getContig() != null && rec.getContig().equals(TARGET_CONTIG) && rec.getAlignmentStart()-250 <= TARGET_START && rec.getAlignmentEnd()+250 >= TARGET_END) {
                    type = "near-target";
                } else {
                    type = "off-target";
                }
            }

            incrementSummaryCounts(type);
            rec.setAttribute("CR",type);
            writer.addAlignment(rec);


            progress.record(rec);
        }

        writer.close();
        it.close();
        CloserUtil.close(reader);

        for(Map.Entry<String, Integer> e : summaryCounts.entrySet()) {
            System.out.println(e.getKey() + " " + e.getValue());
        }
        return 0;
    }

    private void incrementSummaryCounts(String type ) {
        if (summaryCounts.containsKey(type)) {
            summaryCounts.put(type, summaryCounts.get(type) + 1);
        } else {
            summaryCounts.put(type, 1);
        }

    }
}

