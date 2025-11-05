import scipy as sp
import numpy as np

dim_hilbert = 2**10
dim_krylov = 80

def lanczos(A, v, m):
    """
    Algoritmo Lanczos (matrici simmetriche)

    Parametri
    ----------
    A : (n, n) matrice reale simmetrica
        
    v : (n,) vettore di partenza
        
    m : intero, dimensione dello spazio di Krylov
        
    Ritorna
    -------
    V : (n, m) base ortonormale
        
    T : (m, m) matrice tridiagonale che approssima la matrice di partenza
        
    """
    n = A.shape[0]
    V = np.zeros((n, m), dtype=A.dtype)
    T = np.zeros((m, m), dtype=A.dtype)

    v = v / np.linalg.norm(v)
    V[:, 0] = v
    w = A @ v
    alpha = np.vdot(v, w)
    T[0, 0] = alpha
    w = w - alpha * v
    beta = np.linalg.norm(w)

    for j in range(1, m):
        if beta < 1e-12:            #se converge prima si ferma
            return V[:, :j], T[:j, :j]
        v_next = w / beta
        V[:, j] = v_next
        w = A @ v_next - beta * V[:, j - 1]
        alpha = np.vdot(v_next, w)
        w = w - alpha * v_next
        T[j, j] = alpha
        T[j, j - 1] = T[j - 1, j] = beta
        beta = np.linalg.norm(w)

    return V, T


#inizializzo la matrice (simmetrica per velocità e affidabilità) e il vettore iniziale (la funzione lo normalizza)
start_mat = np.random.randint(-3,3, (dim_hilbert, dim_hilbert))
start_mat = (start_mat + start_mat.T)  
start_mat = start_mat.astype(float)
start_vec = np.random.randint(-3,3, size=dim_hilbert)

#print("Starting matrix:\n", start_mat)
#print("Starting vector:\n", start_vec)

krylov_basis, tridiag = lanczos(start_mat, start_vec, dim_krylov)

#print("Orthonormal Krylov basis V:\n", krylov_basis)
#print("\nCheck orthonormality (VᵀV):\n", krylov_basis.T @ krylov_basis)
#print("\nTridiagonal matrix T:\n", tridiag)

#autovalori esatti e approssimati
approx_eigenvals = np.sort(sp.linalg.eigvalsh(tridiag))
exact_eigenvals = np.sort(sp.linalg.eigvalsh(start_mat))
diff = exact_eigenvals[0]-approx_eigenvals[0]
print("\nExact ground state:\n", exact_eigenvals[0])
print("\nApproximate ground state:\n", approx_eigenvals[0])
print("\nError:\n", diff)
print("\nKrylov space dimension:\n", dim_krylov)


	






