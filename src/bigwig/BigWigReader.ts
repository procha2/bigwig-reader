import { DataLoader, BufferedDataLoader, DataMissingError, FileFormatError } from "../loader/DataLoader";
import { BinaryParser } from "../util/BinaryParser";
import { loadHeaderData, HeaderData, FileType } from "./BigWigHeaderReader";
import { loadSequenceRecord, loadSequence, SequenceRecord, streamSequence } from "./TwoBitHeaderReader";
import { inflate } from "pako";
import { Stream, Readable, Writable, Duplex } from "stream";
import { start } from "repl";

export interface BigWigData {
    chr: string,
    start: number,
    end: number,
    value: number
}

export interface BigBedData {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    cdStart?: number,
    cdEnd?: number,
    color?: string,
    exons?: Array<BigBedExon>
}

export interface BigBedDataNarrowPeak {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    signalValue?: number,
    pValue?: number,
    qValue?: number,
    peak?: number,
}

export interface BigBedDataBroadPeak {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    signalValue?: number,
    pValue?: number,
    qValue?: number,
}

export interface BigBedDataRNAElement {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string, // + or - or . for unknown
    level?: number, // Expression level such as RPKM or FPKM. Set to -1 for no data
    signif?: number, // Statistical significance such as IDR. Set to -1 for no data
    score2?: number, // Additional measurement/count e.g. number of reads. Set to 0 for no data
}

export interface BigBedDataMethyl {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    thickStart?: number, // Start of where display should be thick (start codon)
    thickEnd?: number, // End of where display should be thick (stop codon)
    reserved?: number, // Color value R,G,B
    readCount?: number, // Number of reads or coverage
    percentMeth?: number // Percentage of reads that show methylation at this position in the genome
}

export interface BigBedDataTssPeak {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    count?: number, // Count of reads mapping to this peak
    gene_id?: string, // Gene identifier
    gene_name?: string, // Gene name
    tss_id?: string // TSS identifier,
    peak_cov?: string, // base by base read coverage of the peak
}

export interface BigBedDataIdrPeak {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    localIDR?: number, // Local IDR value
    globalIDR?: number, // Global IDR value
    rep1_chromStart?: number, // Start position in chromosome of replicate 1 peak
    rep1_chromEnd?: number, // End position in chromosome of replicate 1 peak
    rep1_count?: number, // Count (used for ranking) replicate 1
    rep2_chromStart?: number, // Start position in chromosome of replicate 2 peak
    rep2_chromEnd?: number, // End position in chromosome of replicate 2 peak
    rep2_count?: number, // Count (used for ranking) replicate 2
}

export interface BigBedDataIdrRankedPeak {
    chr: string,
    start: number,
    end: number,
    name?: string,
    score?: number,
    strand?: string,
    signalValue?: number, // Measurement of enrichment for the region for merged peaks
    pValue?: number, // p-value of merged peak
    qValue?: number, // q-value of merged peak
    summit?: number, // Summit of merged peak
    localIDR?: number, // Local IDR value, which is -log10(local IDR value)
    globalIDR?: number, // Global IDR value, which is -log10(global IDR value)
    chromStart1?: number, // Start position in chromosome of peak 1
    chromEnd1?: number, // End position in chromosome of peak 1
    signalValue1?: number, // Signal measure from peak 1
    summit1?: number, // Summit of peak 1
    chromStart2?: number, // Start position in chromosome of peak 2
    chromEnd2?: number, // End position in chromosome of peak 2
    signalValue2?: number, // Signal measure from peak 2
    summit2?: number, // Summit of peak 2
}

export interface BigBedExon {
    start: number,
    end: number
}

export interface BigZoomData {
    chr: string,
    start: number,
    end: number,
    validCount: number,
    minVal: number,
    maxVal: number,
    sumData: number,
    sumSquares: number
}

interface RPLeafNode {
    startChrom: number;
    startBase: number;
    endChrom: number;
    endBase: number;
    dataOffset: number;
    dataSize: number;
}

const IDX_MAGIC = 0x2468ACE0;
const RPTREE_HEADER_SIZE = 48;
const RPTREE_NODE_LEAF_ITEM_SIZE = 32;
const RPTREE_NODE_CHILD_ITEM_SIZE = 24;
const DEFAULT_BUFFER_SIZE = 512000;

/**
 * Main class for dealing with reading BigWig and BigBed files.
 */
export class BigWigReader {

    private cachedHeader?: HeaderData;
    private cachedSequenceRecords: { [name: string]: SequenceRecord } = {};

    /**
     * @param dataLoader Provided class that deals with fetching data from the file via http, local file, ftp, etc...
     * @param bufferSize Size of the buffer used for fetching data. Used to optimistically read more data than is 
     *      needed for each read of the tree that stores data to avoid round trips. The trade-off is potentially reading 
     *      more data than you need to vs making more round trips.
     */
    constructor(private dataLoader: DataLoader, private bufferSize: number = DEFAULT_BUFFER_SIZE) { }

    /**
     * Gets the type of the underlying file.
     */
    async fileType(): Promise<FileType> {
	    let header: HeaderData = await this.getHeader();
	    return header.fileType;
    }
    
    /**
     * Method for getting all header data for dataLoader's file. Data is loaded on demand and cached for subsequent requests.
     */
    async getHeader(): Promise<HeaderData> {
        if (!this.cachedHeader) {
            this.cachedHeader = await loadHeaderData(this.dataLoader);
        }
        return this.cachedHeader;
    }

    /**
     * Method for getting a sequence record from a 2bit sequence file. This method is not valid for bigWig or bigBed files.
     *
     * @param chrom the name of the chromosome or other sequence to retrieve.
     */
    async getSequenceRecord(chrom: string): Promise<SequenceRecord> {
	    let header: HeaderData = await this.getHeader();
	    if (header.fileType !== FileType.TwoBit) throw new FileFormatError("getSequenceRecord is not valid on " + header.fileType + " files.");
	    if (!this.cachedSequenceRecords[chrom]) {
            this.cachedSequenceRecords[chrom] = await loadSequenceRecord(this.dataLoader, header, chrom);
        }
	    return this.cachedSequenceRecords[chrom];
    }
    
    /**
     * Method for reading unzoomed wig data from BigWig files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     * @param zoomLevelIndex The ZoomLevelHeader.index from the zoom level you want to read from. 
     */
    async readBigWigData(startChrom: string, startBase: number, endChrom: string, 
            endBase: number): Promise<Array<BigWigData>> {
        return this.readData<BigWigData>(startChrom, startBase, endChrom, endBase, 
            (await this.getHeader()).common!.fullIndexOffset, decodeWigData);
    }

    /**
     * Method for streaming unzoomed wig data from BigWig files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     * @param zoomLevelIndex The ZoomLevelHeader.index from the zoom level you want to read from. 
     */
    async streamBigWigData(startChrom: string, startBase: number, endChrom: string, 
            endBase: number): Promise<Readable> {
        return this.streamData<BigWigData>(startChrom, startBase, endChrom, endBase, 
            (await this.getHeader()).common!.fullIndexOffset, decodeWigData);
    }

    /**
     * Method for reading unzoomed bed data from BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedData(startChrom: string, startBase: number, endChrom: string, 
            endBase: number): Promise<Array<BigBedData>> {
        return this.readData<BigBedData>(startChrom, startBase, endChrom, endBase, 
            (await this.getHeader()).common!.fullIndexOffset, decodeBedData);
    }

    /**
     * Method for reading unzoomed bed data from BigBedNarrowPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataNarrowPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataNarrowPeak>> {
    return this.readData<BigBedDataNarrowPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataNarrowPeak);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataBroadPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataBroadPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataBroadPeak>> {
    return this.readData<BigBedDataBroadPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataBroadPeak);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataRNAElement files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataRNAElement(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataRNAElement>> {
    return this.readData<BigBedDataRNAElement>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataRNAElement);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataMethyl files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataMethyl(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataMethyl>> {
    return this.readData<BigBedDataMethyl>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataMethyl);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataTssPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataTssPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataTssPeak>> {
    return this.readData<BigBedDataTssPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataTssPeak);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataIdrPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataIdrPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataIdrPeak>> {
    return this.readData<BigBedDataIdrPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataIdrPeak);
    }

    /**
     * Method for reading unzoomed bed data from BigBedDataIdrRankedPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async readBigBedDataIdrRankedPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Array<BigBedDataIdrRankedPeak>> {
    return this.readData<BigBedDataIdrRankedPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataIdrRankedPeak);
    }

    /**
     * Method for streaming unzoomed bed data from BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedData(startChrom: string, startBase: number, endChrom: string, 
            endBase: number): Promise<Readable> {
        return this.streamData<BigBedData>(startChrom, startBase, endChrom, endBase, 
            (await this.getHeader()).common!.fullIndexOffset, decodeBedData);
    }

    /**
     * Method for streaming unzoomed bed data from BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataBroadPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataBroadPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataBroadPeak);
    }

    /**
     * Method for streaming unzoomed bed data from BigBedDataNarrowPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataNarrowPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataNarrowPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataNarrowPeak);
    }

    /**
     * Method for streaming unzoomed bed data from BigBedDataRNAElement files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataRNAElement(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataRNAElement>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataRNAElement);
    }

    /**
     * Method for streaming unzoomed bed data from BigBedDataMethyl files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataMethyl(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataMethyl>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataMethyl);
    }

    /**
     * Method for streaming unzoomed bed data from BigBedDataTssPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataTssPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataTssPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataTssPeak);
    }

    /**
     * Method for streaming unzoomed bed data from BigBedDataIdrPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataIdrPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataIdrPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataIdrPeak);
    }
    
    /**
     * Method for streaming unzoomed bed data from BigBedDataIdrRankedPeak files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     */
    async streamBigBedDataIdrRankedPeak(startChrom: string, startBase: number, endChrom: string, 
        endBase: number): Promise<Readable> {
    return this.streamData<BigBedDataIdrRankedPeak>(startChrom, startBase, endChrom, endBase, 
        (await this.getHeader()).common!.fullIndexOffset, decodeBigBedDataIdrRankedPeak);
    }

    /**
     * Method for reading Two Bit sequence data from TwoBit files.
     *
     * @param chrom the chromosome from which to read.
     * @param startBase the starting base.
     * @param endBase the ending base.
     */
    async readTwoBitData(chrom: string, startBase: number, endBase: number): Promise<string> {
	    const sequence: SequenceRecord = await this.getSequenceRecord(chrom);
        return loadSequence(this.dataLoader, this.cachedHeader!, sequence, startBase, endBase);
    }

    /**
     * Method for reading Two Bit sequence data from TwoBit files.
     *
     * @param chrom the chromosome from which to read.
     * @param startBase the starting base.
     * @param endBase the ending base.
     */
    async streamTwoBitData(chrom: string, startBase: number, endBase: number, chunkSize: number = 1024): Promise<Readable> {
        const sequence: SequenceRecord = await this.getSequenceRecord(chrom);
        return streamSequence(this.dataLoader, this.cachedHeader!, sequence, startBase, endBase, chunkSize);
    }

    /**
     * Method for reading zoomed data from BigWig and BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     * @param zoomLevelIndex index of the zoom level. You can call getHeader() for a list of these values under HeaderData.zoomLevelHeaders.
     */
    async readZoomData(startChrom: string, startBase: number, endChrom: string, endBase: number, 
            zoomLevelIndex: number): Promise<Array<BigZoomData>> {
        const header = await this.getHeader();
        if (undefined == header.zoomLevelHeaders || !(zoomLevelIndex in header.zoomLevelHeaders)) {
            throw new FileFormatError("Given zoomLevelIndex not found in zoom level headers.");
        }
        const treeOffset = header.zoomLevelHeaders[zoomLevelIndex].indexOffset;
        return this.readData<BigZoomData>(startChrom, startBase, endChrom, endBase, 
            treeOffset, decodeZoomData);
    }

    /**
     * Method for streaming zoomed data from BigWig and BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromose
     * @param endBase Ending base pair
     * @param zoomLevelIndex index of the zoom level. You can call getHeader() for a list of these values under HeaderData.zoomLevelHeaders.
     */
    async streamZoomData(startChrom: string, startBase: number, endChrom: string, endBase: number, 
            zoomLevelIndex: number): Promise<Readable> {
        const header = await this.getHeader();
        if (undefined == header.zoomLevelHeaders || !(zoomLevelIndex in header.zoomLevelHeaders)) {
            throw new FileFormatError("Given zoomLevelIndex not found in zoom level headers.");
        }
        const treeOffset = header.zoomLevelHeaders[zoomLevelIndex].indexOffset;
        return this.streamData<BigZoomData>(startChrom, startBase, endChrom, endBase, 
            treeOffset, decodeZoomData);
    }

    /**
     * Method containing all the shared functionality for reading BigWig and BigBed files.
     * 
     * @param startChrom Starting chromosome
     * @param startBase Starting base pair
     * @param endChrom Ending chromosome
     * @param endBase Ending base pair
     * @param treeOffset Location of the R+ tree that stores the data we're interested.
     * @param decodeFunction 
     */
    private async loadData<T>(startChrom: string, startBase: number, endChrom: string, endBase: number, 
                treeOffset: number, streamMode: boolean, decodeFunction: DecodeFunction<T>, 
                loadFunction: LoadFunction<T>): Promise<void> {
        const header = await this.getHeader();
        if (undefined == header.chromTree) {
            throw new FileFormatError("No chromosome tree found in file header.");
        }
        const startChromIndex: number = header.chromTree.chromToId[startChrom];
        const endChromIndex: number = header.chromTree.chromToId[endChrom];
        if (undefined == startChromIndex) {
            throw new DataMissingError(startChrom);
        }
        if (undefined == endChromIndex) {
            throw new DataMissingError(endChrom);
        }

        // Load all leaf nodes within given chr / base bounds for the R+ tree used for actually storing the data.
        const bufferedLoader = new BufferedDataLoader(this.dataLoader, this.bufferSize, streamMode);
        const magic = new BinaryParser(await bufferedLoader.load(treeOffset, RPTREE_HEADER_SIZE)).getUInt();
        if (IDX_MAGIC !== magic) {
            throw new FileFormatError(`R+ tree not found at offset ${treeOffset}`);
        }
        const rootNodeOffset = treeOffset + RPTREE_HEADER_SIZE;
        const leafNodes: Array<RPLeafNode> = await loadLeafNodesForRPNode(bufferedLoader, header.littleEndian, rootNodeOffset, 
            startChromIndex, startBase, endChromIndex, endBase);

        // Iterate through filtered leaf nodes, load the data, and decode it
        for (const leafNode of leafNodes) {
            let leafData = new Uint8Array(await bufferedLoader.load(leafNode.dataOffset, leafNode.dataSize));
            if (header.common!.uncompressBuffSize > 0) {
                leafData = inflate(leafData);
            }
            let leafDecodedData = decodeFunction(leafData.buffer as ArrayBuffer, startChromIndex, startBase, endChromIndex, 
                endBase, header.chromTree.idToChrom);
            loadFunction(leafDecodedData);
        }
    }

    private async readData<T>(startChrom: string, startBase: number, endChrom: string, endBase: number, 
            treeOffset: number, decodeFunction: DecodeFunction<T>): Promise<Array<T>> {
        const data: Array<T> = [];
        const load: LoadFunction<T> = (d: T[]) => data.push(...d);
        await this.loadData(startChrom, startBase, endChrom, endBase, treeOffset, false, decodeFunction, load);
        return data;
    };

    private async streamData<T>(startChrom: string, startBase: number, endChrom: string, endBase: number, 
            treeOffset: number, decodeFunction: DecodeFunction<T>): Promise<Readable> {
        const stream = new Readable({ objectMode: true, read() {} });
        const load: LoadFunction<T> = (d: T[]) => {
            d.forEach((el) => stream.push(el));
        };
        await this.loadData(startChrom, startBase, endChrom, endBase, treeOffset, true, decodeFunction, load);
        stream.push(null);
        return stream;
    }

}

/**
 * Recursively load a list of R+ tree leaf nodes for the given node (by file offset) within given chr / base bounds.
 * 
 * @param bufferedLoader Buffered data loader used to load the node data.
 * @param rpNodeOffset Offset for the start of the R+ tree node
 * @param startChromIndex starting chromosome index used for filtering
 * @param startBase starting base used for filtering
 * @param endChromIndex ending chromosome index used for filtering
 * @param startBase ending base used for filtering
 * @returns List of simple representations of leaf nodes for the given node offset.
 */
async function loadLeafNodesForRPNode(bufferedLoader: BufferedDataLoader, littleEndian: boolean, rpNodeOffset: number, startChromIndex: number,
        startBase: number, endChromIndex: number, endBase: number): Promise<Array<RPLeafNode>> {
    const nodeHeaderData: ArrayBuffer = await bufferedLoader.load(rpNodeOffset, 4);
    const nodeHeaderParser = new BinaryParser(nodeHeaderData, littleEndian);
    const isLeaf = 1 === nodeHeaderParser.getByte();
    nodeHeaderParser.position++; // Skip reserved space
    const count = nodeHeaderParser.getUShort();

    const nodeDataOffset = rpNodeOffset + 4;
    const bytesRequired = count * (isLeaf ? RPTREE_NODE_LEAF_ITEM_SIZE : RPTREE_NODE_CHILD_ITEM_SIZE);
    const nodeData: ArrayBuffer = await bufferedLoader.load(nodeDataOffset, bytesRequired);

    let leafNodes: Array<RPLeafNode> = [];
    const nodeDataParser = new BinaryParser(nodeData, littleEndian);
    for (let i = 0; i < count; i++) {
        const nodeStartChr = nodeDataParser.getInt();
        const nodeStartBase = nodeDataParser.getInt();
        const nodeEndChr = nodeDataParser.getInt();
        const nodeEndBase = nodeDataParser.getInt();
        // If this node overlaps with the chr / base range provided
        const overlaps: boolean = ((endChromIndex > nodeStartChr) || (endChromIndex == nodeStartChr && endBase >= nodeStartBase)) &&
            ((startChromIndex < nodeEndChr) || (startChromIndex == nodeEndChr && startBase <= nodeEndBase));
        if (isLeaf) {
            const leafNode: RPLeafNode = {
                startChrom: nodeStartChr,
                startBase: nodeStartBase,
                endChrom: nodeEndChr,
                endBase: nodeEndBase,
                dataOffset: nodeDataParser.getLong(),
                dataSize: nodeDataParser.getLong()
            };
            if (overlaps) {
                leafNodes.push(leafNode);
            }
        } else {
            const childOffset = nodeDataParser.getLong();
            if (overlaps) {
                leafNodes.push(... await loadLeafNodesForRPNode(bufferedLoader, littleEndian, childOffset, startChromIndex, startBase, endChromIndex, endBase));
            }
        }
    }

    return leafNodes;
}

type DecodeFunction<T> = (data: ArrayBuffer, startChromIndex: number, startBase: number, endChromIndex: number,
    endBase: number, chromDict: Array<string>) => Array<T>;

type LoadFunction<T> = (data: Array<T>) => void;

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBedData(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
        filterEndBase: number, chromDict: Array<string>): Array<BigBedData> {
    const decodedData: Array<BigBedData> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedData = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseFloat(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.cdStart = parseInt(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.cdEnd = parseInt(tokens[4]);
        }
        if (tokens.length > 5 && tokens[5] !== "." && tokens[5] !== "0") {
            let color: string;
            if (tokens[5].includes(",")) {
                color = tokens[5].startsWith("rgb") ? tokens[5] : "rgb(" + tokens[5] + ")";
            } else {
                color = tokens[5];
            }
            entry.color = color;
        }
        if (tokens.length > 8) {
            const exonCount = parseInt(tokens[6]);
            const exonSizes = tokens[7].split(',');
            const exonStarts = tokens[8].split(',');
            const exons: Array<BigBedExon> = [];

            for (var i = 0; i < exonCount; i++) {
                const eStart = startBase + parseInt(exonStarts[i]);
                const eEnd = eStart + parseInt(exonSizes[i]);
                exons.push({ start: eStart, end: eEnd });
            }

            entry.exons = exons;
        }
        decodedData.push(entry);
    }

    return decodedData;
}


/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataNarrowPeak(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataNarrowPeak> {
const decodedData: Array<BigBedDataNarrowPeak> = [];
const binaryParser = new BinaryParser(data);

const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
while (binaryParser.remLength() >= minSize) {
    const chromIndex = binaryParser.getInt();
    const chrom = chromDict[chromIndex];
    const startBase = binaryParser.getInt();
    const endBase = binaryParser.getInt();
    const rest = binaryParser.getString();

    if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
        continue;
    } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
        break;
    }

    const entry: BigBedDataNarrowPeak = {
        chr: chrom,
        start: startBase,
        end: endBase
    }

    let tokens = rest.split("\t");
    if (tokens.length > 0) {
        entry.name = tokens[0];
    }
    if (tokens.length > 1) {
        entry.score = parseFloat(tokens[1]);
    }
    if (tokens.length > 2) {
        entry.strand = tokens[2];
    }
    if (tokens.length > 3) {
        entry.signalValue = parseInt(tokens[3]);
    }
    if (tokens.length > 4) {
        entry.pValue = parseInt(tokens[4]);
    }
    if (tokens.length > 5) {
        entry.qValue = parseInt(tokens[5]);
    }
    if (tokens.length > 6) {
        entry.peak = parseInt(tokens[6]);
    }

    decodedData.push(entry);
}

return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataBroadPeak(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataBroadPeak> {
    const decodedData: Array<BigBedDataBroadPeak> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataBroadPeak = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseFloat(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.signalValue = parseInt(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.pValue = parseInt(tokens[4]);
        }
        if (tokens.length > 5) {
            entry.qValue = parseInt(tokens[5]);
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataRNAElement(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataRNAElement> {
    const decodedData: Array<BigBedDataRNAElement> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataRNAElement = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseFloat(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.level = parseFloat(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.signif = parseFloat(tokens[4]);
        }
        if (tokens.length > 5) {
            entry.score2 = parseFloat(tokens[5]);
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataMethyl(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataMethyl> {
    const decodedData: Array<BigBedDataMethyl> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataMethyl = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseInt(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.thickStart = parseInt(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.thickEnd = parseInt(tokens[4]);
        }
        if (tokens.length > 5) {
            entry.reserved = parseInt(tokens[5]);
        }
        if (tokens.length > 6) {
            entry.readCount = parseInt(tokens[6]);
        }
        if (tokens.length > 7) {
            entry.percentMeth = parseInt(tokens[7]);
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataTssPeak(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataTssPeak> {
    const decodedData: Array<BigBedDataTssPeak> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataTssPeak = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseFloat(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.count = parseFloat(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.gene_id = tokens[4];
        }
        if (tokens.length > 5) {
            entry.gene_name = tokens[5];
        }
        if (tokens.length > 6) {
            entry.tss_id = tokens[6];
        }
        if (tokens.length > 7) {
            entry.peak_cov = tokens[7];
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataIdrPeak(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataIdrPeak> {
    const decodedData: Array<BigBedDataIdrPeak> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataIdrPeak = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseInt(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.localIDR = parseFloat(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.globalIDR = parseFloat(tokens[4]);
        }
        if (tokens.length > 5) {
            entry.rep1_chromStart = parseInt(tokens[5]);
        }
        if (tokens.length > 6) {
            entry.rep1_chromEnd= parseInt(tokens[6]);
        }
        if (tokens.length > 7) {
            entry.rep1_count = parseFloat(tokens[7]);
        }
        if (tokens.length > 8) {
            entry.rep2_chromStart = parseInt(tokens[8]);
        }
        if (tokens.length > 9) {
            entry.rep2_chromEnd = parseInt(tokens[9]);
        }
        if (tokens.length > 10) {
            entry.rep2_chromEnd = parseFloat(tokens[10]);
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary bed data
 * 
 * @param data Raw bed data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeBigBedDataIdrRankedPeak(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
    filterEndBase: number, chromDict: Array<string>): Array<BigBedDataIdrRankedPeak> {
    const decodedData: Array<BigBedDataIdrRankedPeak> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 3 * 4 + 1;    // Minimum # of bytes required for a bed record
    while (binaryParser.remLength() >= minSize) {
        const chromIndex = binaryParser.getInt();
        const chrom = chromDict[chromIndex];
        const startBase = binaryParser.getInt();
        const endBase = binaryParser.getInt();
        const rest = binaryParser.getString();

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
            break;
        }

        const entry: BigBedDataIdrRankedPeak = {
            chr: chrom,
            start: startBase,
            end: endBase
        }

        let tokens = rest.split("\t");
        if (tokens.length > 0) {
            entry.name = tokens[0];
        }
        if (tokens.length > 1) {
            entry.score = parseInt(tokens[1]);
        }
        if (tokens.length > 2) {
            entry.strand = tokens[2];
        }
        if (tokens.length > 3) {
            entry.signalValue = parseFloat(tokens[3]);
        }
        if (tokens.length > 4) {
            entry.pValue = parseFloat(tokens[4]);
        }
        if (tokens.length > 5) {
            entry.qValue = parseFloat(tokens[5]);
        }
        if (tokens.length > 6) {
            entry.summit= parseInt(tokens[6]);
        }
        if (tokens.length > 7) {
            entry.localIDR = parseFloat(tokens[7]);
        }
        if (tokens.length > 8) {
            entry.globalIDR = parseInt(tokens[8]);
        }
        if (tokens.length > 9) {
            entry.chromStart1 = parseInt(tokens[9]);
        }
        if (tokens.length > 10) {
            entry.chromEnd1 = parseInt(tokens[10]);
        }
        if (tokens.length > 11) {
            entry.signalValue1 = parseFloat(tokens[11]);
        }
        if (tokens.length > 12) {
            entry.summit1 = parseFloat(tokens[12]);
        }
        if (tokens.length > 13) {
            entry.chromStart2 = parseInt(tokens[13]);
        }
        if (tokens.length > 14) {
            entry.chromEnd2 = parseInt(tokens[14]);
        }
        if (tokens.length > 15) {
            entry.signalValue2 = parseFloat(tokens[15]);
        }
        if (tokens.length > 16) {
            entry.summit2 = parseFloat(tokens[16]);
        }

        decodedData.push(entry);
    }

    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary unzoomed wig data
 * 
 * @param data Raw unzoomed wig data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeWigData(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
        filterEndBase: number, chromDict: Array<string>): Array<BigWigData> {
    const decodedData: Array<BigWigData> = [];
    const binaryParser = new BinaryParser(data);

    const chromIndex = binaryParser.getInt();
    const chrom = chromDict[chromIndex];
    let startBase = binaryParser.getInt();
    let endBase = binaryParser.getInt();
    const itemStep = binaryParser.getInt();
    const itemSpan = binaryParser.getInt();
    const type = binaryParser.getByte();
    const reserved = binaryParser.getByte();
    let itemCount = binaryParser.getUShort();

    if (chromIndex < filterStartChromIndex || chromIndex > filterEndChromIndex) {
        return decodedData;
    }

    while (itemCount-- > 0) {
        let value: number;
        if (1 === type) {
            // Data is stored in Bed Graph format
            startBase = binaryParser.getInt();
            endBase = binaryParser.getInt();
            value = binaryParser.getFloat();
        } else if (2 === type) {
            // Data is stored in Variable Step format
            startBase = binaryParser.getInt();
            value = binaryParser.getFloat();
            endBase = startBase + itemSpan;
        } else {
            // Data is stored in Fixed Step format.
            value = binaryParser.getFloat();
            endBase = startBase + itemSpan;
        }

	if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && startBase >= filterEndBase)) {
	    break; // past the end of the range; exit
	} else if (!(chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && endBase < filterStartBase))) {
	    decodedData.push({
		chr: chrom,
		start: startBase,
		end: endBase,
		value: value
            }); // this is within the range (i.e. not before the first requested base); add this datapoint
        }

	if (1 !== type && 2 !== type) {
	    // data is stored in Fixed Step format
	    // only increment the start base once the last entry has been pushed
	    startBase += itemStep;
	}
    }
    return decodedData;
}

/**
 * Extract useful data from sections of raw big binary zoom data
 * 
 * @param data Raw zoomed wig data
 * @param filterStartChromIndex starting chromosome index used for filtering
 * @param filterStartBase starting base used for filtering
 * @param filterEndChromIndex ending chromosome index used for filtering
 * @param filterEndBase ending base used for filtering
 * @param chromDict dictionary of indices used by the file to chromosome names, conveniently stored as an array.
 */
function decodeZoomData(data: ArrayBuffer, filterStartChromIndex: number, filterStartBase: number, filterEndChromIndex: number,
        filterEndBase: number, chromDict: Array<string>): Array<BigZoomData> {
    const decodedData: Array<BigZoomData> = [];
    const binaryParser = new BinaryParser(data);

    const minSize = 8 * 4;   // Minimum # of bytes required for a zoom record
    while (binaryParser.remLength() > minSize) {
        const chromIndex = binaryParser.getInt();
        const decodedZoomData: BigZoomData = {
            chr: chromDict[chromIndex],
            start: binaryParser.getInt(),
            end: binaryParser.getInt(),
            validCount: binaryParser.getInt(),
            minVal: binaryParser.getFloat(),
            maxVal: binaryParser.getFloat(),
            sumData: binaryParser.getFloat(),
            sumSquares: binaryParser.getFloat()
        };

        if (chromIndex < filterStartChromIndex || (chromIndex === filterStartChromIndex && decodedZoomData.end < filterStartBase)) {
            continue;
        } else if (chromIndex > filterEndChromIndex || (chromIndex === filterEndChromIndex && decodedZoomData.start >= filterEndBase)) {
            break;
        }
        decodedData.push(decodedZoomData);
    }
    return decodedData;
}
