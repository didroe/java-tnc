/*
 * See the "LICENSE" file for the full license governing this code.
 */
package didroe.tnc;

/**
 * Criteria that caused the minimisation to complete
 *
 * @author Did
 */
public enum CompletionReason {

    /**
     * Local minima reached (|pg| ~= 0)
     */
    TNC_LOCALMINIMUM("Local minima reached (|pg| ~= 0)"),
    /**
     * Converged (|f_n-f_(n-1)| ~= 0)
     */
    TNC_FCONVERGED("Converged (|f_n-f_(n-1)| ~= 0)"),
    /**
     * Converged (|x_n-x_(n-1)| ~= 0)
     */
    TNC_XCONVERGED("Converged (|x_n-x_(n-1)| ~= 0)");

    private final String message;

    CompletionReason(String message) {
        this.message = message;
    }

    public String getMessage() {
        return message;
    }
}
