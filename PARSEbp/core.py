import os, math, time, argparse, re, tempfile, subprocess, sys, stat
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np

class ScoreResult:
    def __init__(self, score_arr, decoylist, rejected_pdblist):
        self.scores = np.array(score_arr, dtype=float)
        self.decoylist = list(decoylist)
        self.rejected_pdblist = list(rejected_pdblist)

    def save(self, filename):
        # Determine output path
        if os.path.isdir(filename):
            outputfilename = os.path.join(filename, "score.txt")
        else:
            outputfilename = filename

        with open(outputfilename, "w") as f:
            # Write valid decoys and scores
            for name, score in zip(self.decoylist, self.scores):
                f.write(f"{name}\t{score:.3f}\n")

            # Write rejected decoys with placeholder score
            for rej_path in self.rejected_pdblist:
                rej_name = os.path.basename(rej_path)
                f.write(f"{rej_name}\t-\n")

    def getScore(self, pdbname):
        """Return the score for the given decoy name."""
        if pdbname in self.rejected_pdblist:
            print(f"{pdbname} is not scored")
            return 0
        try:
            idx = self.decoylist.index(pdbname)
            return round(float(self.scores[idx]), 3)
        except ValueError:
            raise KeyError(f"Decoy '{pdbname}' not found in decoy list.")

    def top1(self):
        """
        Return the top scored decoy(s) and the score.
        If multiple decoys share the same top score, return all.
        Returns:
            (list[str], float)
        """
        max_score = np.max(self.scores)
        top_decoys = [name for name, s in zip(self.decoylist, self.scores) if s == max_score]
        return top_decoys, round(float(max_score), 3)

    def topN(self, N=5):
        """
        Return a dictionary of the top N scoring decoys sorted by score descending.
        Args:
            N (int): number of top entries to return
        Returns:
            dict[str, float]
        """
        if N <= 0:
            raise ValueError("N must be a positive integer.")

        sorted_indices = np.argsort(self.scores)[::-1]  # descending
        top_indices = sorted_indices[:N]
        return {self.decoylist[i]: round(float(self.scores[i]), 3) for i in top_indices}


class parsebp:
    def __init__(self):
        self.mode = 1
        self.num_threads = 50
        self.base_dir = os.path.dirname(os.path.abspath(__file__))

        self.target_sequence = ""
        self.pdb_dir=""
        self.length = 0
        
        self.decoylist = []
        self.pdblist = []
        self.rejected_pdblist = []

        self.usalign = self._get_executable("USalign")
        self.mcann = self._get_executable("MC-Annotate") 

    def _get_executable(self, name: str) -> str:
        """Get absolute path to bundled executable and ensure +x permissions."""
        exe_path = os.path.join(self.base_dir, "bin", name)

        if not os.path.exists(exe_path):
            raise FileNotFoundError(f"Executable not found: {exe_path}")

        # Ensure executable permission
        st = os.stat(exe_path)
        os.chmod(exe_path, st.st_mode | stat.S_IEXEC)

        return exe_path

    def set_target_sequnece(self, target_sequence):

        self.target_sequence = target_sequence

    def set_mode(self, mode):

        self.mode = mode

    def set_parallel_threads(self, num_threads):

        self.num_threads = num_threads

    def load_pdbs(self, dir_path):
        """
        Load all PDBs from a directory.
        If self.target_sequence == "", all PDBs are accepted.
        Otherwise, only those with matching sequence are kept.
        """

        self.pdb_dir = dir_path
        all_pdb_files = [f for f in os.listdir(dir_path) if f.endswith(".pdb")]

        seq_lengths = []

        self.decoylist = []
        self.pdblist = []
        self.rejected_pdblist = []

        for fname in tqdm(all_pdb_files, desc="Loading PDBs"):
            full_path = os.path.join(dir_path, fname)

            seq, length = self.get_sequence_and_length(full_path)
            seq_lengths.append(length)

            if self.target_sequence == "" or seq == self.target_sequence:
                self.decoylist.append(fname)
                self.pdblist.append(full_path)
            else:
                self.rejected_pdblist.append(fname)

        # Determine L
        if self.target_sequence == "":
            self.length = max(seq_lengths) if seq_lengths else 0
        else:
            self.length = len(self.target_sequence)

        if len(self.pdblist) == 0:
            sys.exit("No PDBs to score in specified input directory")
        if len(self.pdblist) == 1:
            sys.exit("You must provide atleast 2 PDBs to score in the specified input directory")
        
        print(f"Loaded {len(self.pdblist)} matching PDBs to score, rejected {len(self.rejected_pdblist)}")
        
        if self.rejected_pdblist:
            print("Rejected pdbs:")
            for f in self.rejected_pdblist:
                print(f)

        
        
    def get_sequence_and_length(self, pdbfilename):
        """
        Extract the sequence and its length from a PDB file.
        Returns (sequence_string, length).
        """
        sequence = []
        seen_residues = set()

        with open(pdbfilename, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("TER"):
                    break
                if not line.startswith("ATOM"):
                    continue

                res_name = line[17:20].strip()    # residue name (e.g., A, G, C, U)
                res_id = int(line[22:26].strip()) # residue number

                if res_id not in seen_residues:
                    seen_residues.add(res_id)
                    sequence.append(res_name)

        seq_str = "".join(sequence)
        return seq_str, len(seq_str)

    def parseTMscore(self, filename):
        """Parse USalign log file and return the main TM-score."""
        tmscore, tmcount = None, -1
        with open(filename, 'r', encoding='UTF-8') as file:
            for line in file:
                if "TM-score" in line:
                    tmcount += 1
                    if tmcount % 3 == 0:
                        tokenlist = line.split("=")
                        tokenlist2 = tokenlist[1].split()
                        tmscore = float(tokenlist2[0])
                        break
        return tmscore


    def calcTM(self, args):
        """Run USalign and return TM-score."""
        nativefilename, fullpdbfilename, counter = args
        with tempfile.TemporaryDirectory() as tmpdir:
            logfilename = os.path.join(tmpdir, f"tmlog{counter}.txt")
            cmd = [f"{self.usalign}", fullpdbfilename, nativefilename]
            with open(logfilename, "w") as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.DEVNULL, check=False)
            return self.parseTMscore(logfilename)


    def calcINF(self, args):
        pred_map, ref_map, _ = args
        mask = np.triu(np.ones_like(pred_map, dtype=bool), k=1)
        pred_1 = pred_map[mask] == 1
        nat_1 = ref_map[mask] == 1

        tp = np.sum(pred_1 & nat_1)
        fp = np.sum(pred_1 & ~nat_1)
        fn = np.sum(~pred_1 & nat_1)

        if tp + fp == 0 or tp + fn == 0:
            return 0.0
        return math.sqrt((tp / (tp + fp)) * (tp / (tp + fn)))


    def compute_pairwise_matrix_upper(self, modelslist, func, num_threads=50):
        """
        Compute NxN pairwise matrix (upper triangle only), using threads for I/O/subprocess tasks.
        """
        numdecoys = len(modelslist)
        matrix = np.zeros((numdecoys, numdecoys))
        
        # Build argument list and corresponding index pairs
        args_list = [(modelslist[i], modelslist[j], f"{i}_{j}")
                    for i in range(numdecoys) for j in range(i+1, numdecoys)]
        pair_index = [(i, j)
                    for i in range(numdecoys) for j in range(i+1, numdecoys)]

        scores = []
        total_batches = math.ceil(len(args_list) / num_threads)
        total_tasks = len(args_list)

        if func.__name__ == "calcTM":
            desc_str = "Computing pairwise TM-scores"
        else:
            desc_str = "Computing pairwise INF scores"
        
        with tqdm(total=total_tasks, desc=desc_str, dynamic_ncols=True) as pbar:
            for b in range(total_batches):
                batch = args_list[b*num_threads:(b+1)*num_threads]

                batch_scores = [None] * len(batch)
                with ThreadPoolExecutor(max_workers=min(num_threads, len(batch))) as executor:
                    future_to_idx = {executor.submit(func, arg): idx for idx, arg in enumerate(batch)}
                    for future in as_completed(future_to_idx):
                        idx = future_to_idx[future]
                        try:
                            batch_scores[idx] = future.result()
                        except Exception as e:
                            batch_scores[idx] = 0.0
                            print(f"[Warning] Thread failed: {e}")
                        pbar.update(1)

                scores.extend(batch_scores)

        for k, (i, j) in enumerate(pair_index):
            matrix[i, j] = matrix[j, i] = scores[k]

        np.fill_diagonal(matrix, 1.0)
        
        return matrix


    def getEdgeFromBasePair(self, args):
        fullmodelpath, L, counter = args

        # Check chain presence from PDB (first ATOM line)
        chain_present = 0
        with open(fullmodelpath, 'r', encoding='UTF-8') as ifile:
            for line in ifile:
                if line.startswith("ATOM"):
                    chain = line[21:22]
                    if chain != " ":
                        chain_present = 1
                    break

        # Run MC-Annotate and capture output directly (no file written!)
        output = subprocess.check_output(
            [f"{self.mcann}", fullmodelpath],
            universal_newlines=True
        )

        ssmap = np.zeros((L, L), dtype=int)

        counter=0

        # Parse the output lines
        lines = output.splitlines()
        for idx, line in enumerate(lines):
            if line.startswith("Base-pairs"):
                for line in lines[idx+1:]:
                    tokenlist = line.split(":")
                    if not tokenlist:  # skip bad lines
                        continue
                    tstr = tokenlist[0]

                    hiphen = tstr.find("-")
                    all_ind = [m.start() for m in re.finditer("'", tstr)]

                    if len(all_ind) == 0:
                        node1 = int(tstr[chain_present:hiphen]) - 1
                        node2 = int(tstr[hiphen+2-(1-chain_present):-1]) - 1
                    else:
                        node1 = int(tstr[all_ind[1]+1:hiphen]) - 1
                        node2 = int(tstr[all_ind[3]+1:-1]) - 1

                    ssmap[node1][node2] = 1
                    ssmap[node2][node1] = 1
                    counter += 1

        return ssmap
        

    def extract2D(self, pdblist, L, num_threads=50):
        args_list = [(pdblist[i], L, i) for i in range(len(pdblist))]
        ss_maps = [None] * len(args_list)

        total_batches = math.ceil(len(args_list) / num_threads)
        total_tasks = len(args_list)

        with tqdm(total=total_tasks, desc="Extracting base pairs", dynamic_ncols=True) as pbar:
            for batch_index in range(total_batches):
                batch = args_list[batch_index * num_threads: (batch_index + 1) * num_threads]

                with ThreadPoolExecutor(max_workers=min(num_threads, len(batch))) as executor:
                    futures = {executor.submit(self.getEdgeFromBasePair, arg): arg[2] for arg in batch}

                    for f in as_completed(futures):
                        idx = futures[f]
                        try:
                            ss_maps[idx] = f.result()
                        except Exception as e:
                            ss_maps[idx] = None
                            print(f"[Warning] base pair extraction failed at index {idx}: {e}")
                        pbar.update(1)

        return ss_maps
    
    def score(self):
        
        start_t = time.time()

        # # Step 1: build pairwise matrix
        matrix = self.compute_pairwise_matrix_upper(self.pdblist, func = self.calcTM, num_threads = self.num_threads)

        if self.mode:

            ssmaps = self.extract2D(self.pdblist, self.length, num_threads = self.num_threads)
        
            matrixSS = self.compute_pairwise_matrix_upper(ssmaps, func = self.calcINF, num_threads = self.num_threads)
        
            matrix = matrix * matrixSS

        scores = (matrix.sum(axis=1) - np.diag(matrix)) / (matrix.shape[0] - 1)

        time_elp = time.time() - start_t

        print(f"Scored {len(self.pdblist)} decoys\tTime taken = {time_elp:.2f} seconds")

        return ScoreResult(scores, self.decoylist, self.rejected_pdblist)
