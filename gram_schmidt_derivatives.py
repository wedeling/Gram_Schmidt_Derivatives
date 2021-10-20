"""
Class for computing the derivatives of Gram-Schmidt vectors.

Author: Wouter Edeling

"""

import numpy as np

class Gram_Schmidt_Derivatives:
    """
    Class for computing the derivatives of Gram-Schmidt vectors
    """

    def __init__(self, Q):
        """
        Inintialize Gram_Schmidt_Derivative object. This will already
        compute the Gram-Schmidt orthogonalization of Q, yielding a list
        w = [w_1, ..., w_d] of orthogonal vectors.

        Parameters
        ----------
        Q : array size (D, d)
            An array of d independent (non-orthogonal) vectors of size D.

        Returns
        -------
        None.

        """
        self.D = Q.shape[0]
        self.d = Q.shape[1]
        # list for orthogonal vectors w_i
        self.w = []
        # list for orthonormal vectors w_i/||w_i||_2
        self.normed_w = []
        # list for norms ||w_i||_2
        self.norm_w = []
        # store q vectos in list as well
        self.q = [Q[:, i].reshape([self.D, 1]) for i in range(self.d)]
        # compute w_i and w_i/||w_i||_2
        self.gram_schmidt()
        # compute derivatives of w_i and w_i/||w_i||_2
        self.derivatives()

    def set_Q(self, Q):
        """
        Overwrite the independent non-orthogonal basis vectors q_i and recompute
        Gram-Schmidt.

        Parameters
        ----------
        Q : array size (D, d)
            An array of d independent (non-orthogonal) vectors of size D.

        Returns
        -------
        None.

        """
        self.__init__(Q)

    def get_orthogonal_vectors(self):
        """
        Returns the Gram-Schmidt vectors w_i.

        Returns
        -------
        W : array, size = (D,d)
            Array containing the d orthogonal vectors w_i.

        """
        W = np.zeros([self.D, self.d])
        for i in range(self.d):
            W[:, i] = self.w[i].flatten()
        return W

    def get_orthonormal_vectors(self):
        """
        Returns the Gram-Schmidt vectors w_i / ||w_i||_2.

        Returns
        -------
        W : array, size = (D,d)
            Array containing the d orthogonal vectors w_i.

        """
        W = np.zeros([self.D, self.d])
        for i in range(self.d):
            W[:, i] = self.normed_w[i].flatten()
        return W

    def get_w_derivatives(self):
        """
        Get all the Gram-Schmidt derivatives dw_i/ dq_k.

        Returns
        -------
        dw_dq : dict
            The Gram-Schmidt derivatives. The keys of the dictionary are (i, k)
            indices corresponding to dw_i/ dq_k, with i,k = 1, ..., d.

        """
        dw_dq = {}
        #loop over dw_i, i = 1, 2,...,d
        for i in range(self.d):
            #loop over all dq_j, j = 1,...,d
            for j in range(self.d):
                if i >= j:
                    #convert 2D index (i,j) to scalar index idx
                    idx = np.ravel_multi_index([i, j], dims=(self.d, self.d))
                    dw_dq[(i+1, j+1)] = self.grads[idx]
        return dw_dq

    def get_normed_w_derivatives(self):
        """
        Get all the Gram-Schmidt derivatives d(w_i/||w_i||_2)/ dq_k.

        Returns
        -------
        dw_dq : dict
            The Gram-Schmidt derivatives. The keys of the dictionary are (i, k)
            indices corresponding to d(w_i/||w_i||_2)/ dq_k, with i,k = 1, ..., d.

        """
        dw_dq = {}
        #loop over dw_i, i = 1, 2,...,d
        for i in range(self.d):
            #loop over all dq_j, j = 1,...,d
            for j in range(self.d):
                if i >= j:
                    #convert 2D index (i,j) to scalar index idx
                    idx = np.ravel_multi_index([i, j], dims=(self.d, self.d))
                    dw_dq[(i+1, j+1)] = self.grads_normed[idx]
        return dw_dq

    def gram_schmidt(self):
        """
        Gram-Schmidt orthogonalization. Non-orthogonal basis vectors are stored internally.
        To overwrite these use set_Q() method.

        Returns
        -------
        None.

        """

        print("Computing Gram-Schmidt orthogonalization...")
        # first orthogonal vector is just q_0
        w_0 = np.copy(self.q[0])
        self.w.append(w_0)
        # also compute normalized w
        self.norm_w.append(np.linalg.norm(w_0))
        self.normed_w.append(w_0 / self.norm_w[0])
        # orthogonalize the remaining d-1 vectors
        for i in range(1, self.d):
            q_i = np.copy(self.q[i])
            # start with q_i
            current_w = q_i
            for j in range(1, i + 1):
                w_j = self.w[j-1]
                # subtract orthogonal projection of q_i on previous w_j vectors
                current_w -= np.dot(w_j.T, q_i)/np.dot(w_j.T, w_j) * w_j
            # store w, its norm and the normed vector w / ||w||_2
            self.w.append(current_w)
            self.norm_w.append(np.linalg.norm(current_w))
            self.normed_w.append(current_w / self.norm_w[i])
        print("done.")

    def compute_D_ij(self, w_j, q_i):
        """
        Compute the D_ij matrices which appear in the gradients of the
        Gram-Schmidt vectors w_i wrt the original unnormalized vectors q_j.
        These matrices are only defined for i unequal to j.

        D_ij = 1/(w_j^Tw_j) [w_jq_i^T - (2w_j^Tq_i)/(w_j^Tw_j)*w_jw_j^T + w_j^Tq_i*I]
        where I is the D x D identity matrix.

        Parameters
        ----------
        w_j : array, shape (D,1)
            The orthogonal, but not yet normalized, vector w_i, i=1,...,d.
        q_i : array, shape (D,1)
            An orginal, non-orthonogal, vector q_j, j=1,...,d.

        Returns
        -------
        D_ij : array
            The D x D matrix described above.

        """

        return (np.dot(w_j, q_i.T) -
                (2.0*np.dot(w_j.T, q_i) / np.dot(w_j.T, w_j)) * np.dot(w_j, w_j.T) +
                np.dot(w_j.T, q_i) * np.eye(self.D)) / np.dot(w_j.T, w_j)

    def derivatives(self):
        """
        Compute the derivatives of the normalized Gram-Schmidt vectors w_i with respect to
        the original non-orthogobal basis vectors q_k. Also computes the derivatives
        of the normed vectors w_i / ||w_i||_2.

        Returns
        -------
        grads : array, size (d**2, D, D)
            The d**2 derivatives dw_i / dq_k for i,k = 1, ... , d.
        grads_normed : array, size (d**2, D, D)
            The d**2 derivatives of the normed Gram-Schmidt vector:
                d(w_i / ||w_i||_2) / dq_k for i,k = 1, ... , d..

        """
        print("Computing Gram-Schmidt derivatives...")

        d = self.d
        D = self.D
        w = self.w
        q = self.q
        norm_w = self.norm_w

        #gradients
        grads = np.zeros([d**2, D, D])
        #normed gradients
        grads_normed = np.zeros([d**2, D, D])
        #D_ij matrices
        D_ij = np.zeros([d**2, D, D])
        #for w_1, the gradient and D_ij matrices have a simple form. Compute these
        #outside loop.
        I_D = np.eye(D)
        D_ij[0] = I_D
        grads[0] = I_D
        grads_normed[0] = I_D/norm_w[0] - np.dot(w[0], w[0].T)/norm_w[0]**3

        #the D_ij matrices must be precomputed
        #loop over dw_i, i = 2,...,d
        for i in range(1, d):
            #loop over all dq_j, j = 1,...,d
            for j in range(d):
                #convert 2D index (i,j) to scalar index idx
                idx = np.ravel_multi_index([i, j], dims=(d,d))
                if i >= j:
                    D_ij[idx] = self.compute_D_ij(w[j], q[i])

        #construct the d^2 gradient matrices dw_i / dq_k
        for i in range(1, d):
            #the matrix by which the gradient d_wi / dq_k must be premultiplied
            norm_mat = I_D/norm_w[i] - np.dot(w[i], w[i].T)/norm_w[i]**3
            for k in range(d):
                #convert 2D index (i,j) to scalar index idx
                idx = np.ravel_multi_index([i, k], dims=(d,d))
                #due to GS, the gradient dw_i / dq_k is zero for k > i, so compute if i >= k
                if i >= k:
                    #a 'normal' gradient dw_i / dq_i
                    if i == k:
                        # index of dw_{i-1} / dq_{i-1}
                        idx1 = np.ravel_multi_index([i-1, i-1], dims = (d,d))
                        grad_ik = grads[idx1] - np.dot(w[i-1], w[i-1].T)/np.dot(w[i-1].T, w[i-1])
                    #a 'shear' gradient dw_i / dq_k where i is not k
                    else:
                        grad_ik = 0.0
                        #loop over j=1,...,i - 1
                        for j in range(i):
                            #index of D_{ij}
                            idx2 = np.ravel_multi_index([i, j], dims = (d,d))
                            #index of dw_j / dq_k
                            idx3 = np.ravel_multi_index([j, k], dims = (d,d))
                            grad_ik -= np.dot(D_ij[idx2], grads[idx3])
                    #store dw_i / dq_k
                    grads[idx] = grad_ik
                    #compute the derivates of w_i / ||w_i||_2
                    grads_normed[idx] = np.dot(norm_mat, grad_ik)

        self.grads = grads
        self.grads_normed = grads_normed
        print('done.')
