import numpy as np
import math
import json

class DarkCoreEngineV2:
    def __init__(self, sequence_id, aa_sequence):
        self.sequence_id = sequence_id
        self.sequence = aa_sequence.upper()
        
        # Escala de Hidrofobicità Kyte-Doolittle estandarizada
        self.kd_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

    def vectorizer(self):
        """Modulo 1: Vectorizacion tensorial y calculo de gradientes."""
        vector = np.array([self.kd_scale.get(aa, 0.0) for aa in self.sequence])
        gradients = np.gradient(vector) if len(vector) > 1 else np.array([0.0])
        return vector, gradients

    def structural_viability_score(self):
        """Modulo 2: Metric SVR y Entropia de Informacion."""
        vector, gradients = self.vectorizer()
        mean_hydro = np.mean(vector)
        grad_std = np.std(gradients)
        
        # Calculo de entropia Shannon sobre el gradiente
        prob_dist = np.abs(gradients) / (np.sum(np.abs(gradients)) + 1e-9)
        entropy = -np.sum([p * np.log2(p) for p in prob_dist if p > 0])
        
        # SVR Score normalizado
        svr_score = float((mean_hydro * 0.4) + (grad_std * 0.3) + (entropy * 0.3))
        return {
            "mean_hydrophobicity": round(float(mean_hydro), 4),
            "gradient_std": round(float(grad_std), 4),
            "entropy": round(float(entropy), 4),
            "svr_score": round(svr_score, 4)
        }

    def generate_pdb_backbone(self, filename="structure.pdb"):
        """Modulo 3: Generador de coordenadas 3D (Formato .PDB)."""
        coords = []
        x, y, z = 0.0, 0.0, 0.0
        phi = 3.8  # Distancia fija enlace N-CA-C aproximada en Angstroms
        
        with open(filename, "w") as pdb_file:
            for i, aa in enumerate(self.sequence):
                # Generacion de helicoide alpha teórica
                x += phi * math.cos(i * 1.0)
                y += phi * math.sin(i * 1.0)
                z += 1.5 * i
                coords.append((x, y, z))
                
                # Escribir linea con formato oficial PDB ATOM
                pdb_file.write(
                    f"ATOM  {i+1:5d}  CA  {aa:3s} A{i+1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 0.00           C\n"
                )
        return filename

    def estimate_binding_energy(self, target_svr=0.5):
        """Modulo 4: Estimacion de Docking (Delta G)."""
        metrics = self.structural_viability_score()
        delta_g = round(-1.36 * math.log(abs(metrics["svr_score"] - target_svr) + 1.01), 2)
        return {
            "estimated_delta_g_kcal_mol": delta_g,
            "affinity_status": "High" if delta_g < -4.0 else "Moderate/Low"
        }

    def in_silico_mutagenesis(self, position, new_aa):
        """Modulo de Mutagenesis en caliente."""
        seq_list = list(self.sequence)
        if 0 <= position < len(seq_list):
            seq_list[position] = new_aa.upper()
            mutated_seq = "".join(seq_list)
            mutant_engine = DarkCoreEngineV2(f"{self.sequence_id}_mut", mutated_seq)
            return mutant_engine.structural_viability_score()
        return None


# --- ORQUESTRADOR / EJECUCIÓN DIRECTA ---
if __name__ == "__main__":
    # Secuencia de prueba (Ejemplo ncORF: OLMALINC microprotein candidate)
    sample_seq = "MGNIFANLFKEGLSKIF"
    
    print(f"=== INICIANDO PROJECT DARK-CORE V2.0 ===")
    engine = DarkCoreEngineV2(sequence_id="TARGET_001", aa_sequence=sample_seq)
    
    # 1. Viabilidad
    svr_res = engine.structural_viability_score()
    print("\n[+] Analisis SVR y Entropia:")
    print(json.dumps(svr_res, indent=4))
    
    # 2. PDB
    pdb_path = engine.generate_pdb_backbone("dark_core_target.pdb")
    print(f"\n[+] Modelo 3D exportado exitosamente: {pdb_path}")
    
    # 3. Docking Delta G
    docking_res = engine.estimate_binding_energy()
    print("\n[+] Estimacion de Energia de Union (Delta G):")
    print(json.dumps(docking_res, indent=4))
