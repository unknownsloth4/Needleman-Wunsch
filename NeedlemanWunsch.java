public class NeedlemanWunsch {
	protected final String seq1, seq2;
	public int scoring;
	protected Node[][] matrix;
	
	protected final int matchScore, mismatchScore, gapPenalty;
	
	String atas = "";	// sequence asam amino 1
	String tengah = "";	// Memasangkan sequence asam amino atas(1) & bawah(2)
						// (dg garis kalau sama atau dg blank space kalau beda)
	String bawah = "";	// sequence asam amino 2
	
	public NeedlemanWunsch (String s1, String s2, Penilaian s) {
		seq1 = "\u25A0" + s1;
		seq2 = "\u25A0" + s2;
		
		// menentukan jenis penilaian, apakah asam amino sama, berbeda, atau ada gap
		matchScore = s.matchScore;
		mismatchScore = s.mismatchScore;
		gapPenalty = s.gapPenalty;
		
		// inisiasi baris dan kolom pertama dari matriks dengan nilai gap penalty
		matrix = new Node[seq1.length()][seq2.length()];
		
		for (int i = 0; i < seq1.length(); i++) {
			matrix[i][0] = new Node(i * gapPenalty, i, 0);
		}
		
		for (int i = 0; i < seq2.length(); i++) {
			matrix[0][i] = new Node(i * gapPenalty, 0, i);
		}
	}
	
	// menentukan penilaian apa yg digunakan (miss atau match)
	protected int nilai(int i, int j) {
		if (seq1.charAt(i) == seq2.charAt(j)) {
			return matchScore;
		} else {
			return mismatchScore;
		}
	}
	
	// scoring dg algoritma NW
	// mencari nilai yang terbesar untuk diisikan ke matrix
	protected Node match (int i, int j) {
		Node s1, s2, s3;
		s1 = new Node (matrix[i-1][j-1].nilai + nilai(i, j), i, j);	// diagonal
		s2 = new Node (matrix[i-1][j].nilai + gapPenalty, i, j);	// atas
		s3 = new Node (matrix[i][j-1].nilai + gapPenalty, i, j);	// kiri
		
		Node terbesar =  new Node(Math.max(s1.nilai, Math.max(s2.nilai, s3.nilai)), i, j);
		if (s1.compareTo(terbesar) == 0) {
			terbesar.prev = matrix[i-1][j-1];
		} else if (s2.compareTo(terbesar) == 0) {
			terbesar.prev = matrix[i-1][j];
		}
        else
            terbesar.prev = matrix[i][j-1];

        return terbesar;
	}
	
	// menerapkan algoritma NW ke setiap node pd matrix
	// isi matrix dg nilai yg didapat
	public Node isiMatrix() {
        for (int i = 1; i < seq1.length(); i++) {
            for (int j = 1; j < seq2.length(); j++){
                matrix[i][j] = match(i, j);
            }
        }
        
        // return nilai node terakhir pada matrix utk traceback
        scoring = matrix[seq1.length()-1][seq2.length()-1].nilai;
        return matrix[seq1.length()-1][seq2.length()-1];
    }
    
    // menentukan alignment pada sequence asam amino
    protected Node alignment(Node posisi) {
        while (posisi.prev != null) {

            if (posisi.i - posisi.prev.i == 1 && posisi.j - posisi.prev.j == 1){    // Diagonal
                boolean x = seq1.charAt(posisi.i) == seq2.charAt(posisi.j) ? true : false;
                if(x) {
                	atas = seq1.charAt(posisi.i) + " " + atas;
                    tengah = "|" + " " + tengah;
                    bawah = seq2.charAt(posisi.j) + " " + bawah;
                } else {
                    atas = seq1.charAt(posisi.i) + " " + atas;
                    tengah = " " + " " + tengah;
                    bawah = seq2.charAt(posisi.j) + " " + bawah;
                }
            } else if (posisi.i - posisi.prev.i == 1) {	// Atas
                atas = seq1.charAt(posisi.i) + " " + atas;
                tengah = " " + " " + tengah;
                bawah = "-" + " " + bawah;
            } else if (posisi.j - posisi.prev.j == 1) {	// Kiri
                atas = "-" + " " + atas;
                tengah = " " + " " + tengah;
                bawah = seq2.charAt(posisi.j) + " " + bawah;
            }

            posisi = posisi.prev;
        }

        return posisi;
    }
    
    // traceback dari node terakhir matrix ke awal
    public void traceback() {
        Node posisi = matrix[seq1.length()-1][seq2.length()-1];
        posisi = alignment(posisi);

		// jika path berakhir pada indel (dua sequence asam amino memiliki awal yang berbeda/hang)
		// agar sequence asam amino dan path bisa ditulis lengkap
		// akan kembali ke 0,0
        while (posisi.i != 0 || posisi.j != 0) {
            if (posisi.i != 0 && posisi.j == 0){
                posisi.prev = matrix[posisi.i-1][posisi.j];
                posisi = alignment(posisi);
            }else if (posisi.i == 0 && posisi.j != 0) {
                posisi.prev = matrix[posisi.i][posisi.j-1];
                posisi = alignment(posisi);
            }
        }
		
		// print alignment
        System.out.println(atas);
        System.out.println(tengah);
        System.out.println(bawah);
    }
    
    // mencetak matrix dg nilai akhir yg didapat
    public void printMatrix() {
        System.out.printf("%4s", "\u25A0");
        for (int i = 0; i < matrix[0].length; i++) {
            System.out.printf("%4s", seq2.charAt(i));
        }
        
        System.out.println();
        
        for (int i = 0; i < matrix.length; i++) {
            System.out.printf("%4s", seq1.charAt(i));
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%4d", matrix[i][j].nilai);
            }
            
            System.out.println();
        }
    }
    
    public static void main(String[] args) {
        String seq1 = "ATCG";
        String seq2 = "TCG";

        Penilaian s = new Penilaian(4, 0, -6);             
        NeedlemanWunsch seq = new NeedlemanWunsch(seq1, seq2, s);     

        seq.isiMatrix();
        seq.traceback();
        System.out.print("\n");
        seq.printMatrix();
        
        System.out.println("\nScore : " + seq.scoring + "\n");
    }
}

// inisiasi class Node utk membuat linked list pada matrix    
class Node implements Comparable<Node>{
    int i, j;
    int nilai;
    Node prev;

    public Node(int nilai, int x, int y) {
        this.i = x;
        this.j = y;
        this.nilai = nilai;
        this.prev = null;
    }

    public int compareTo(Node n) {
        return this.nilai - n.nilai;
    }
}

// insiasi parameter utk penilaian alignment
class Penilaian {
    int matchScore, mismatchScore, gapPenalty;

    public Penilaian(int ms, int mms, int gp) {
        matchScore = ms;
        mismatchScore = mms;
        gapPenalty = gp;
    }
}