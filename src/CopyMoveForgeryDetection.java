import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.imageio.ImageIO;

public class CopyMoveForgeryDetection {
    private static class FileImageData {
        public final File file;
        public final BufferedImage image;

        public FileImageData(File file, BufferedImage image) {
            this.file = file;
            this.image = image;
        }
    }

    // Loading the images
    public static List<FileImageData> loadImages (String path){
        List<FileImageData> images = new ArrayList<>();
        File folder = new File(path);
        for (File file : folder.listFiles()) {
            try {
                BufferedImage image = ImageIO.read(file);
                FileImageData data = new FileImageData(file, image);
                images.add(data);
            } catch (Exception e) {
                // handle exception
            }
        }
        return images;
    }

    // Constants for block size and threshold
    public static final int BLOCK_SIZE = 16;
    public static final double THRESHOLD = 10;

    public static void main(String[] args) {
        // Load the CoMoFoD database
        String dir = System.getProperty("user.dir");
        List<FileImageData> images = loadImages(dir + "\\CoMoFoD\\Forged");

        for (FileImageData imageData : images) {
            String currentFile = imageData.file.getName();
            currentFile = currentFile.substring(0, currentFile.lastIndexOf("."));
            BufferedImage image = imageData.image;

            // Create a new file for logging the correlations
            File logFile = new File(dir+"\\Correlations\\log" + currentFile + ".txt");
            BufferedWriter bw = null;

            try {
                bw = new BufferedWriter(new FileWriter(logFile));
                bw.write("Current BLOCK_SIZE = " + BLOCK_SIZE + ", current THRESHOLD = " + THRESHOLD + "\n");
                // Divide the image into blocks
                int[][] blocks = divideIntoBlocks(image);

                // Compute the DCT of each block
                double[][][] dcts = computeDCTs(blocks);

                // Compare DCT coefficients of each block to all other blocks
                for (int i = 0; i < dcts.length; i++) {
                    for (int j = 0; j < dcts[i].length; j++) {
                        // Skip blocks that have already been matched
                        if (blocks[i][j] == -1) {
                            continue;
                        }
                        for (int k = i; k < dcts.length; k++) {
                            for (int l = (k == i) ? j + 1 : 0; l < dcts[k].length; l++) {
                                // Skip blocks that have already been matched
                                if (blocks[k][l] == -1) {
                                    continue;
                                }
                                // Compare DCT coefficients
                                double[] diff = new double[BLOCK_SIZE * BLOCK_SIZE];
                                for (int m = 0; m < BLOCK_SIZE; m++) {
                                    for (int n = 0; n < BLOCK_SIZE; n++) {
                                        diff[m * BLOCK_SIZE + n] = Math.abs(dcts[i][j][m * BLOCK_SIZE + n] - dcts[k][l][m * BLOCK_SIZE + n]);

                                    }
                                }
                                double meanDiff = mean(diff);
                                if (meanDiff < THRESHOLD) {
                                    // Match found
                                    System.out.println(currentFile + " - Match found between block (" + i + ", " + j + ") and block (" + k + ", " + l + ")");
                                    bw.write(currentFile + " - Match found between block (" + i + ", " + j + ") and block (" + k + ", " + l + ")\n");
                                    blocks[i][j] = -1;
                                    blocks[k][l] = -1;
                                }
                            }
                        }
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            } finally {
                try {
                    if (bw != null) {
                        bw.close();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    // Divide an image into blocks
    public static int[][] divideIntoBlocks (BufferedImage image){
        int[][] blocks = new int[image.getWidth() / BLOCK_SIZE][image.getHeight() / BLOCK_SIZE];
        for (int i = 0; i < blocks.length; i++) {
            for (int j = 0; j < blocks[i].length; j++) {
                blocks[i][j] = image.getRGB(i * BLOCK_SIZE, j * BLOCK_SIZE);
            }
        }
        return blocks;
    }

    // Compute the DCT of each block
    public static double[][][] computeDCTs (int[][] blocks){
        double[][][] dcts = new double[blocks.length][blocks[0].length][BLOCK_SIZE * BLOCK_SIZE];
        for (int i = 0; i < blocks.length; i++) {
            for (int j = 0; j < blocks[i].length; j++) {
                // Convert block to grayscale
                double[][] grayBlock = new double[BLOCK_SIZE][BLOCK_SIZE];
                for (int k = 0; k < BLOCK_SIZE; k++) {
                    for (int l = 0; l < BLOCK_SIZE; l++) {
                        int rgb = blocks[i][j];
                        int r = (rgb >> 16) & 0xff;
                        int g = (rgb >> 8) & 0xff;
                        int b = rgb & 0xff;
                        grayBlock[k][l] = (r + g + b) / 3.0;
                    }
                }
                // Compute DCT of grayscale block
                double[][] dct = computeDCT(grayBlock);
                for (int k = 0; k < BLOCK_SIZE; k++) {
                    for (int l = 0; l < BLOCK_SIZE; l++) {
                        dcts[i][j][k * BLOCK_SIZE + l] = dct[k][l];
                    }
                }
            }
        }
        return dcts;
    }

    // Compute the DCT of a block
    public static double[][] computeDCT (double[][] block){
        double[][] dct = new double[BLOCK_SIZE][BLOCK_SIZE];
        for (int u = 0; u < BLOCK_SIZE; u++) {
            for (int v = 0; v < BLOCK_SIZE; v++) {
                double sum = 0.0;
                for (int i = 0; i < BLOCK_SIZE; i++) {
                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        sum += block[i][j] * Math.cos(((2 * i + 1) / (2.0 * BLOCK_SIZE)) * u * Math.PI) * Math.cos(((2 * j + 1) / (2.0 * BLOCK_SIZE)) * v * Math.PI);
                    }
                }
                dct[u][v] = sum;
            }
        }
        return dct;
    }

    // Calculate mean of array
    public static double mean (double[] arr){
        double sum = 0.0;
        for (double val : arr) {
            sum += val;
        }
        return sum / arr.length;
    }
}