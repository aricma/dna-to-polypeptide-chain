#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include <map>


enum Nucleotide {
    Adenine,
    Thymine,
    Uracil,
    Cytosine,
    Guanine,
};


enum Protein {
    Alanine,
    Valine,
    Arginine,
    Serine,
    Lysine,
    Asparagine,
    Threonine,
    Methionine,
    Isoleucine,
    Glutamine,
    Histidine,
    Proline,
    Leucine,
    Tryptophan,
    StopCodon,
    Cysteine,
    Tyrosine,
    Phenylalanine,
    Glycine,
    GlutamicAcid,
    AsparticAcid,
};


Nucleotide dna_character_to_nucleotide(char dna_character) {
    if (dna_character == 'A') {
        return Adenine;
    } else if (dna_character == 'T') {
        return Thymine;
    } else if (dna_character == 'U') {
        return Uracil;
    } else if (dna_character == 'C') {
        return Cytosine;
    } else if (dna_character == 'G') {
        return Guanine;
    } else {
        std::string dna_character_as_string;
        dna_character_as_string = dna_character;
        throw std::invalid_argument(
            dna_character_as_string + " <- is not a nucleotide character! Use A, T, U, C or G instead."
        );
    }
}


char nucleotide_to_dna_character(Nucleotide nucleotide) {
    if (nucleotide == Adenine) return 'A';
    if (nucleotide == Thymine) return 'T';
    if (nucleotide == Uracil) return 'U';
    if (nucleotide == Cytosine) return 'C';
    // if (nucleotide == Guanine)
    else return 'G';
}

std::string protein_to_string(Protein protein) {
    if (protein == Alanine) return "Alanine";
    if (protein == Valine) return "Valine";
    if (protein == Arginine) return "Arginine";
    if (protein == Serine) return "Serine";
    if (protein == Lysine) return "Lysine";
    if (protein == Asparagine) return "Asparagine";
    if (protein == Threonine) return "Threonine";
    if (protein == Methionine) return "Methionine";
    if (protein == Isoleucine) return "Isoleucine";
    if (protein == Glutamine) return "Glutamine";
    if (protein == Histidine) return "Histidine";
    if (protein == Proline) return "Proline";
    if (protein == Leucine) return "Leucine";
    if (protein == Tryptophan) return "Tryptophan";
    if (protein == StopCodon) return "StopCodon";
    if (protein == Cysteine) return "Cysteine";
    if (protein == Tyrosine) return "Tyrosine";
    if (protein == Phenylalanine) return "Phenylalanine";
    if (protein == Glycine) return "Glycine";
    if (protein == GlutamicAcid) return "GlutamicAcid";
    if (protein == AsparticAcid) return "AsparticAcid";
    else {
        throw std::invalid_argument( "protein_to_string does only support the Protein enum" );
    }
}


std::string get_nucleotide_tuple(Nucleotide nucleotides [3]) {
    std::string nucleotide_tuple;

    for (int i = 0; i < 3; i++) {
        std::string nucleotide_as_string;
        nucleotide_as_string = nucleotide_to_dna_character(nucleotides[i]);
        nucleotide_tuple = nucleotide_tuple + nucleotide_as_string;
    }

    return nucleotide_tuple;
}


Protein nucleotides_to_protein(Nucleotide nucleotides [3]) {
    std::string nucleotide_tuple = get_nucleotide_tuple(nucleotides);

    if (nucleotide_tuple == "AGA") return Arginine;
    if (nucleotide_tuple == "AGG") return Arginine;

    if (nucleotide_tuple == "AGC") return Serine;
    if (nucleotide_tuple == "AGU") return Serine;

    if (nucleotide_tuple == "AAA") return Lysine;
    if (nucleotide_tuple == "AAG") return Lysine;

    if (nucleotide_tuple == "AAC") return Asparagine;
    if (nucleotide_tuple == "AAU") return Asparagine;

    if (nucleotide_tuple == "ACG") return Threonine;
    if (nucleotide_tuple == "ACA") return Threonine;
    if (nucleotide_tuple == "ACC") return Threonine;
    if (nucleotide_tuple == "ACU") return Threonine;

    if (nucleotide_tuple == "AUG") return Methionine;

    if (nucleotide_tuple == "AUA") return Isoleucine;
    if (nucleotide_tuple == "AUC") return Isoleucine;
    if (nucleotide_tuple == "AUU") return Isoleucine;

    if (nucleotide_tuple == "CGA") return Arginine;
    if (nucleotide_tuple == "CGU") return Arginine;
    if (nucleotide_tuple == "CGC") return Arginine;
    if (nucleotide_tuple == "CGG") return Arginine;

    if (nucleotide_tuple == "CAG") return Glutamine;
    if (nucleotide_tuple == "CAA") return Glutamine;

    if (nucleotide_tuple == "CAC") return Histidine;
    if (nucleotide_tuple == "CAU") return Histidine;

    if (nucleotide_tuple == "CCC") return Proline;
    if (nucleotide_tuple == "CCG") return Proline;
    if (nucleotide_tuple == "CCA") return Proline;
    if (nucleotide_tuple == "CCU") return Proline;

    if (nucleotide_tuple == "CUC") return Leucine;
    if (nucleotide_tuple == "CUG") return Leucine;
    if (nucleotide_tuple == "CUA") return Leucine;
    if (nucleotide_tuple == "CUU") return Leucine;

    if (nucleotide_tuple == "UGG") return Tryptophan;

    if (nucleotide_tuple == "UGA") return StopCodon;

    if (nucleotide_tuple == "UGC") return Cysteine;
    if (nucleotide_tuple == "UGU") return Cysteine;

    if (nucleotide_tuple == "UAG") return StopCodon;
    if (nucleotide_tuple == "UAA") return StopCodon;

    if (nucleotide_tuple == "UAC") return Tyrosine;
    if (nucleotide_tuple == "UAU") return Tyrosine;

    if (nucleotide_tuple == "UCA") return Serine;
    if (nucleotide_tuple == "UCU") return Serine;
    if (nucleotide_tuple == "UCC") return Serine;
    if (nucleotide_tuple == "UCG") return Serine;

    if (nucleotide_tuple == "UUA") return Leucine;
    if (nucleotide_tuple == "UUU") return Leucine;

    if (nucleotide_tuple == "UUC") return Phenylalanine;
    if (nucleotide_tuple == "UUG") return Phenylalanine;

    if (nucleotide_tuple == "GAA") return GlutamicAcid;
    if (nucleotide_tuple == "GAG") return GlutamicAcid;

    if (nucleotide_tuple == "GAU") return AsparticAcid;
    if (nucleotide_tuple == "GAC") return AsparticAcid;

    if (nucleotide_tuple == "GUA") return Valine;
    if (nucleotide_tuple == "GUU") return Valine;
    if (nucleotide_tuple == "GUC") return Valine;
    if (nucleotide_tuple == "GUG") return Valine;

    if (nucleotide_tuple == "GCA") return Alanine;
    if (nucleotide_tuple == "GCU") return Alanine;
    if (nucleotide_tuple == "GCC") return Alanine;
    if (nucleotide_tuple == "GCG") return Alanine;

    if (nucleotide_tuple == "GGA") return Glycine;
    if (nucleotide_tuple == "GGU") return Glycine;
    if (nucleotide_tuple == "GGG") return Glycine;
    if (nucleotide_tuple == "GGC") return Glycine;

    else {
        throw std::invalid_argument(
                nucleotide_tuple + " <- is not a nucleotide_tuple or has no matching Protein!"
        );
    }
}


Nucleotide get_complementary_nucleotide(Nucleotide nucleotide) {
    if (nucleotide == Adenine) return Thymine;
    if (nucleotide == Thymine) return Adenine;
    if (nucleotide == Uracil) return Adenine;
    if (nucleotide == Cytosine) return Guanine;
    // if (nucleotide == Guanine)
    else return Cytosine;
}


std::vector<Nucleotide> dna_to_complementary_rna(std::vector<Nucleotide> dna) {
    std::vector<Nucleotide> complementary_rna;

    for(std::vector<int>::size_type i = 0; i != dna.size(); i++) {
        Nucleotide complementary_nucleotide;
        complementary_nucleotide = get_complementary_nucleotide(dna[i]);
        if (complementary_nucleotide == Thymine) {
            complementary_nucleotide = Uracil;
        }
        complementary_rna.push_back(complementary_nucleotide);
    }

    return complementary_rna;
}


std::vector<Nucleotide> read_dna_from_input() {
    std::vector<Nucleotide> dna;

    char dna_character;
    while (std::cin >> dna_character) {
        Nucleotide nucleotide = dna_character_to_nucleotide(dna_character);
        dna.push_back(nucleotide);
    }

    return dna;
}

std::vector<Protein> get_polypeptide_chain(std::vector<Nucleotide> rna) {
    std::vector<Protein> polypeptide_chain;
    Nucleotide next_nucleotide_tuple [3];
    int nucleotide_tuple_counter = 0;

    for(std::vector<int>::size_type i = 0; i != rna.size(); i++) {
        next_nucleotide_tuple[nucleotide_tuple_counter] = rna[i];
        nucleotide_tuple_counter++;
        if (nucleotide_tuple_counter == 3) {
            Protein protein = nucleotides_to_protein(next_nucleotide_tuple);
            polypeptide_chain.push_back(protein);
            nucleotide_tuple_counter = 0;
        }
    }

    return polypeptide_chain;
}


void write_line(std::string line) {
    std::cout << line << std::endl;
}


void write_lines(std::vector<std::string> lines) {
    for(std::vector<int>::size_type i = 0; i != lines.size(); i++) {
        write_line(lines[i]);
    }
}


void write_nucleotides(std::vector<Nucleotide> nucleotides) {
    for(std::vector<int>::size_type i = 0; i != nucleotides.size(); i++) {
        std::string line;
        line = nucleotide_to_dna_character(nucleotides[i]);
        write_line(line);
    }
}

void write_polypeptide_chain(std::vector<Protein> proteins) {
    for(std::vector<int>::size_type i = 0; i != proteins.size(); i++) {
        std::string line;
        line = protein_to_string(proteins[i]);
        write_line(line);
    }
}


int main() {
    std::vector<Nucleotide> dna = read_dna_from_input();
    std::vector<Nucleotide> complementary_rna = dna_to_complementary_rna(dna);
    std::vector<Protein> polypeptide_chain = get_polypeptide_chain(complementary_rna);
    write_polypeptide_chain(polypeptide_chain);
    return 0;
}