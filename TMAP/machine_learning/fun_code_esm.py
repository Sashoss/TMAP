


from Bio.SeqUtils import ProtParam

import numpy as np
import tensorflow as tf



def convert_embedding_to_sequence(esm2_embeddings):
    """Converts esm2 embeddings to an amino acid sequence.

    Args:
        esm2_embeddings: A list of esm2 embeddings, one for each amino acid in the sequence.

    Returns:
        A string representing the amino acid sequence.
    # Example usage
    esm2_embeddings = [0.1, 0.2, 0.3, 0.4, 0.5]
    sequence = convert_embedding_to_sequence(esm2_embeddings)
    print(sequence)
    """
    # Calculate the properties of the amino acid residues
    residues = []
    for embedding in esm2_embeddings:
        residue = ProtParam.ProteinAnalysis(embedding)
        residues.append(residue)

    # Map the residues back to their corresponding amino acid letters
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    residue_to_aa = {residues[i]: amino_acids[i] for i in range(len(amino_acids))}
    sequence = ''.join([residue_to_aa[residue] for residue in residues])

    return sequence





def evolve():

    

    # Load the esm2 embeddings for the input and evolved amino acid sequences
    input_seq_embeddings = np.load("input_seq_esm2_embeddings.npy")
    evolved_seq_embeddings = np.load("evolved_seq_esm2_embeddings.npy")

    # Define the model architecture
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(64, input_shape=(input_seq_embeddings.shape[1],), activation='relu'))
    model.add(tf.keras.layers.Dense(64, activation='relu'))
    model.add(tf.keras.layers.Dense(evolved_seq_embeddings.shape[1], activation='softmax'))

    # Compile the model
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

    # Train the model on the input and evolved amino acid sequence embeddings
    model.fit(input_seq_embeddings, evolved_seq_embeddings, epochs=10, batch_size=32)

    # Use the trained model to predict the evolved amino acid sequence for a new input sequence
    new_input_seq_embedding = np.load("new_input_seq_esm2_embedding.npy")
    predicted_evolved_seq_embedding = model.predict(new_input_seq_embedding)

    # Convert the predicted evolved sequence embedding back into a sequence of amino acids
    predicted_evolved_seq = convert_embedding_to_sequence(predicted_evolved_seq_embedding)

    print(predicted_evolved_seq)
