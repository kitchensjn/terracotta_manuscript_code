import numpy as np
from scipy import linalg
import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct



def _calc_current_pos(id, messages, parents, coal_rate):
    """Calculates current node position as product of child messages

    Parameters
    ----------
    id : int
        ID of node
    messages : np.array
        Messages being passed in tree
    parents : np.array
        Parent IDs for each node
    
    Returns
    -------
    current_pos : np.array
        Probability distribution of node's current position given subtree below
    """

    return np.multiply(coal_rate, np.prod(messages[np.where(parents==id)[0]], axis=0))

def _calc_branch_message(
        current_pos,
        branch_above,
        transition_matrices,
        direction="backward"
    ):
    """Calculates the message to be passed along a branch above specified node

    Parameters
    ----------
    id : int
        ID of node
    current_pos : np.array
        Probability distribution of node's current position given subtree below
    branch_above : np.array
        Branch lengths above each node split across epochs. Shape is #epochs x #nodes.
    direction : string


    Returns
    -------
    current_pos : np.array
        Probability distribution for location of lineage given subtree below. Length is #demes.
    """


    included_epochs = np.where(branch_above > 0)[0]
    for epoch in included_epochs:
        trans_prob = linalg.expm(transition_matrices[epoch]*branch_above[epoch])
        if direction == "backward":
            current_pos = np.matmul(trans_prob, current_pos)
        else:
            current_pos = np.matmul(current_pos, trans_prob)
    return current_pos

def likelihood_of_tree(
        parents,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    ):
    """

    Parameters
    ----------
    parents : np.array
        Parent IDs for each node
    branch_above : np.array
        Branch lengths above each node split across epochs
    node_epoch
    ids_asc_time : np.array
        Nodes IDs in time ascending order
    sample_locations_array : np.array
        Array
    sample_ids : np.array
        Defines order of sample node IDs for sample_locations_array
    transition_matrices : np.array
        Rate matrices for each epoch
    coal_rate

    Returns
    -------
    loglikelihood : float
        Log-likelihood of tree
    """

    num_demes = len(sample_locations_array[0])
    messages = np.ones((len(parents), num_demes), dtype="float64")
    loglikelihood = 0
    for id in ids_asc_time: 
        if id in sample_ids:
            current_pos = sample_locations_array[np.where(sample_ids==id)[0][0]]
        else:
            current_pos = _calc_current_pos(
                id,
                messages,
                parents,
                coal_rates[node_epoch[id]]
            )
        parent = parents[id]
        if parent != -1:
            messages[id] = _calc_branch_message(
                current_pos,
                branch_above[:,id],
                transition_matrices
            )
        else:   # collect roots here
            messages[id] = current_pos
            loglikelihood += np.log(np.sum(current_pos))
    return loglikelihood, messages



ids_asc_time = np.array([0, 1, 2, 3, 4])
parents = np.array([3, 3, 4, 4, -1])
node_times = np.array([0, 0, 0, 2, 4])
branch_above = np.array([[2, 2, 4, 2, 0]])
node_epoch = np.array([0, 0, 0, 0, 0])
sample_ids = np.array([0, 1, 2])
sample_locations_array = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
])

m = 0.01
s = np.array([[0.25, 0.5, 1.0]])
coal_rates = 1/s
transition_matrices = np.array([
    [
        [-(m*(s[0][1]/s[0][0]) + m*(s[0][2]/s[0][0])), m*(s[0][0]/s[0][1]), m*(s[0][0]/s[0][2])],
        [m*(s[0][1]/s[0][0]), -(m*(s[0][0]/s[0][1]) + m*(s[0][2]/s[0][1])), m*(s[0][1]/s[0][2])],
        [m*(s[0][2]/s[0][0]), m*(s[0][2]/s[0][1]), -(m*(s[0][0]/s[0][2]) + m*(s[0][1]/s[0][2]))]
    ]
])

from0_1b = _calc_branch_message(
    np.array([1, 0, 0]),
    np.array([2]),
    transition_matrices
)
from1_1b = _calc_branch_message(
    np.array([0, 1, 0]),
    np.array([2]),
    transition_matrices
)
from2_2b = _calc_branch_message(
    np.array([0, 0, 1]),
    np.array([4]),
    transition_matrices
)

final = 0
for j in range(3):
    total = 0
    for i in range(3):
        node3 = np.array([0, 0, 0])
        node3[i] = 1
        backward = _calc_branch_message(
            node3,
            np.array([2]),
            transition_matrices
        )
        total += from0_1b[i]*from1_1b[i]*backward[j]*from2_2b[j]*coal_rates[0][i]*coal_rates[0][j]
    final += total

loglikelihood, messages = likelihood_of_tree(
    parents,
    branch_above,
    node_epoch,
    ids_asc_time,
    sample_locations_array,
    sample_ids,
    transition_matrices,
    coal_rates
)


for i in range(3):
    total = 0
    for j in range(3):
        node4 = np.array([0, 0, 0])
        node4[j] = 1
        forward = _calc_branch_message(
            node4,
            np.array([1]),
            transition_matrices,
            direction="forward"
        )
        total += from0_1b[i]*from1_1b[i]*forward[i]*from2_2b[j]*coal_rates[0][i]*coal_rates[0][j]


def calc_all_messages(
        parents,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    ):
    """"""

    num_demes = len(sample_locations_array[0])
    messages = np.ones((len(parents)*2, num_demes), dtype="float64")
    for id in ids_asc_time: 
        if id in sample_ids:
            current_pos = sample_locations_array[np.where(sample_ids==id)[0][0]]
        else:
            current_pos = _calc_current_pos(
                id,
                messages,
                parents,
                coal_rates[node_epoch[id]]
            )
        parent = parents[id]
        if parent != -1:
            messages[id] = _calc_branch_message(
                current_pos,
                branch_above[:,id],
                transition_matrices
            )
        else:   # collect roots here
            messages[id] = current_pos
    for id in ids_asc_time[::-1]:
        parent_of = np.where(parents==id)[0]
        for c in range(len(parent_of)):
            alt_children = np.delete(parent_of, c)
            current_pos = np.multiply(coal_rates[node_epoch[id]], np.prod(np.concatenate((messages[[id+len(parents)]], messages[alt_children])), axis=0))
            messages[parent_of[c]+len(parents)] = _calc_branch_message(
                current_pos,
                branch_above[:,parent_of[c]],
                transition_matrices,
                direction="forward"
            )
    return messages

def trace_ancestors(start, parents):
    """
    """

    lineage = []
    a = start
    while a != -1:
        lineage.append(a)
        a = parents[a]
    return np.array(lineage)

def track_lineage_over_time(
        sample,
        times,
        parents,
        node_times,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    ):

    ancestors = trace_ancestors(sample, parents)
    ancestor_times = node_times[ancestors]

    pc_combos = []
    for t in times:
        for i,v in enumerate(ancestor_times):
            if v > t:
                child = ancestors[i-1]
                parent = ancestors[i]
                break
            elif v == t:
                child = ancestors[i]
                parent = ancestors[i]
                break
        pc_combos.append((child, parent))

    messages = calc_all_messages(
        parents,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    )

    positions = np.zeros((len(pc_combos), len(sample_locations_array[0])))
    for element, node_combo in enumerate(pc_combos):
        if node_combo[0] == node_combo[1]:
            if node_combo[0] in sample_ids:
                node_pos = sample_locations_array[np.where(sample_ids==node_combo[0])[0][0]]
            else:
                combined = np.multiply(coal_rates[node_epoch[node_combo[0]]], np.prod(np.concatenate((messages[[node_combo[0]+len(parents)]], messages[np.where(parents==node_combo[0])[0]])), axis=0))
                node_pos = combined / sum(combined)
        else:
            if node_combo[0] in sample_ids:
                child_pos = sample_locations_array[np.where(sample_ids==node_combo[0])[0][0]]
            else:
                incoming_child_messages = messages[np.where(parents==node_combo[0])[0]]
                if len(incoming_child_messages) > 0:
                    combined = np.multiply(coal_rates[node_epoch[node_combo[0]]], np.prod(incoming_child_messages, axis=0))
                    child_pos = combined / sum(combined)
                else:
                    child_pos = np.ones((1,len(sample_locations_array[0])))[0]
            if node_combo[1] in sample_ids:
                parent_pos = sample_locations_array[np.where(sample_ids==node_combo[1])[0][0]]
            else:
                backward_messages = np.where(parents==node_combo[1])[0]
                backward_messages = backward_messages[backward_messages != node_combo[0]]
                incoming_parent_messages = np.concatenate((messages[[node_combo[1]+len(parents)]], messages[backward_messages]))
                if len(incoming_parent_messages) > 0:
                    combined = np.multiply(coal_rates[node_epoch[node_combo[1]]], np.prod(incoming_parent_messages, axis=0))
                    parent_pos = combined / sum(combined)
                else:
                    parent_pos = np.ones((1,len(sample_locations_array[0])))[0]
            
            branch_length_to_child = int(times[element] - node_times[node_combo[0]])
            bl_child = branch_above[:, node_combo[0]].copy()
            bl_parent = bl_child.copy()
            whats_left = branch_length_to_child
            for e in range(len(branch_above)):
                current = bl_parent[e].copy()
                if current >= whats_left:
                    bl_parent[e] -= whats_left
                else:
                    bl_parent[e] = 0
                whats_left -= current
            bl_child -= bl_parent

            outgoing_child_message = _calc_branch_message(
                child_pos,
                bl_child,
                transition_matrices,
                direction="backward"
            )
            outgoing_parent_message = _calc_branch_message(
                parent_pos,
                bl_parent,
                transition_matrices,
                direction="forward"
            )

            node_pos = np.multiply(outgoing_child_message, outgoing_parent_message)
        positions[element] = node_pos / sum(node_pos)

    return positions

sample = 0
times = [0, 1, 2, 3, 4]

track_lineage_over_time(
    sample,
    times,
    parents,
    node_times,
    branch_above,
    node_epoch,
    ids_asc_time,
    sample_locations_array,
    sample_ids,
    transition_matrices,
    coal_rates
)











#print(messages)
#print(np.multiply(coal_rates[0], np.prod(messages[[0, 1, 8]], axis=0)))


exit()

def _calc_all_messages(
        parents,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    ):

    num_demes = len(sample_locations_array[0])
    messages = np.ones((len(parents)*2, num_demes), dtype="float64")
    for id in ids_asc_time:
        # Within this loop, repeatedly implement Equation 2 from the manuscript

        if id in sample_ids:
            current_pos = sample_locations_array[np.where(sample_ids==id)[0][0]]
        else:
            current_pos = np.prod(messages[np.where(parents==id)[0]], axis=0)
            #current_pos = np.multiply(coal_rates[node_epoch[id]], np.prod(messages[np.where(parents==id)[0]], axis=0))

        #current_pos = current_pos / np.sum(current_pos)
        
        messages[id] = _calc_branch_message(
            current_pos,
            branch_above[:,id],
            transition_matrices,
            direction="backward"
        )

    for id in ids_asc_time[::-1]:
        parent_of = np.where(parents==id)[0]
        for c in range(len(parent_of)):
            alt_children = np.delete(parent_of, c)
            
            current_pos = np.prod(np.concatenate((messages[[id+len(parents)]], messages[alt_children])), axis=0)
            #current_pos = np.multiply(coal_rates[node_epoch[id]], np.prod(np.concatenate((messages[[id+len(parents)]], messages[alt_children])), axis=0))

            #current_pos = current_pos / np.sum(current_pos)

            messages[parent_of[c]+len(parents)] = _calc_branch_message(
                current_pos,
                branch_above[:,parent_of[c]],
                transition_matrices,
                direction="forward"
            )

    return messages



messages = _calc_all_messages(
    parents,
    branch_above,
    node_epoch,
    ids_asc_time,
    sample_locations_array,
    sample_ids,
    transition_matrices,
    coal_rates
)

print(messages)
exit()

print(np.multiply(np.prod(np.concatenate((messages[[4+len(parents)]], messages[np.where(parents==4)[0]])), axis=0), coal_rates))


exit()






exit()




def _calc(
        parents,
        branch_above,
        node_epoch,
        ids_asc_time,
        sample_locations_array,
        sample_ids,
        transition_matrices,
        coal_rates
    ):

    num_demes = len(sample_locations_array[0])




messages = _calc_all_messages(
    parents,
    branch_above,
    node_epoch,
    ids_asc_time,
    sample_locations_array,
    sample_ids,
    transition_matrices,
    coal_rates
)


exit()











ids_asc_time = np.array([0, 1, 2, 3, 4, 5, 6])
parents = np.array([4, 4, 5, 6, 5, 6, -1])
branch_above = np.array([[1, 1, 1, 1, 1, 1, 0]])
sample_locations_array = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]
])
sample_ids = np.array([0, 1, 2, 3])

transition_matrices = np.array([
    [
        [-3, 2, 1],
        [1, -2, 1],
        [1, 1, -2]
    ]
])