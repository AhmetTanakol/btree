#ifndef INCLUDE_MODERNDBS_BTREE_H
#define INCLUDE_MODERNDBS_BTREE_H

#include <cstddef>
#include <cstring>
#include <functional>
#include <optional>
#include "moderndbs/buffer_manager.h"
#include "moderndbs/defer.h"
#include "moderndbs/segment.h"

namespace moderndbs {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct BTree : public Segment {
    struct Node {
        /// The level in the tree.
        uint16_t level;

        /// The number of children.
        uint16_t count;

        // Constructor
        Node(uint16_t level, uint16_t count)
            : level(level), count(count) {}

        /// Is the node a leaf node?
        bool is_leaf() const { return level == 0; }
    };

    struct InnerNode: public Node {
        /// The capacity of a node.
        /// TODO think about the capacity that the nodes have.
        static constexpr uint32_t kCapacity = 42;

        /// The keys.
        KeyT keys[kCapacity];

        /// The children.
        uint64_t children[kCapacity];

        /// Constructor.
        InnerNode() : Node(0, 0) {}

        /// Get the index of the first key that is not less than than a provided key.
        /// @param[in] key          The key that should be searched.
        std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
            // HINT
        }

        /// Insert a key.
        /// @param[in] key          The separator that should be inserted.
        /// @param[in] split_page   The id of the split page that should be inserted.
        void insert(const KeyT &key, uint64_t split_page) {
            // HINT
        }

        /// Split the node.
        /// @param[in] buffer       The buffer for the new page.
        /// @return                 The separator key.
        KeyT split(std::byte* buffer) {
            // HINT
        }

        /// Returns the keys.
        /// Can be implemented inefficiently as it's only used in the tests.
        std::vector<KeyT> get_key_vector() {
            // TODO
        }

        /// Returns the child page ids.
        /// Can be implemented inefficiently as it's only used in the tests.
        std::vector<uint64_t> get_child_vector() {
            // TODO
        }
    };

    struct LeafNode: public Node {
        /// The capacity of a node.
        /// TODO think about the capacity that the nodes have.
        static constexpr uint32_t kCapacity = 42;

        /// The keys.
        KeyT keys[kCapacity];

        /// The values.
        ValueT values[kCapacity];

        /// Constructor.
        LeafNode() : Node(0, 0) {}

        /// Insert a key.
        /// @param[in] key          The key that should be inserted.
        /// @param[in] value        The value that should be inserted.
        void insert(const KeyT &key, const ValueT &value) {
            // HINT
        }

        /// Erase a key.
        void erase(const KeyT &key) {
            // HINT
        }

        /// Split the node.
        /// @param[in] buffer       The buffer for the new page.
        /// @return                 The separator key.
        KeyT split(std::byte* buffer) {
            // HINT
        }

        /// Returns the keys.
        /// Can be implemented inefficiently as it's only used in the tests.
        std::vector<KeyT> get_key_vector() {
            // TODO
        }

        /// Returns the values.
        /// Can be implemented inefficiently as it's only used in the tests.
        std::vector<ValueT> get_value_vector() {
            // TODO
        }
    };

    /// The root.
    std::optional<uint64_t> root;

    /// Next page id.
    /// You don't need to worry about about the page allocation.
    /// (Neither fragmentation, nor persisting free-space bitmaps)
    /// Just increment the next_page_id whenever you need a new page.
    uint64_t next_page_id;

    /// Constructor.
    BTree(uint16_t segment_id, BufferManager &buffer_manager)
        : Segment(segment_id, buffer_manager) {
        // TODO
    }

    /// Lookup an entry in the tree.
    /// @param[in] key      The key that should be searched.
    std::optional<ValueT> lookup(const KeyT &key) {
        // TODO
    }

    /// Erase an entry in the tree.
    /// @param[in] key      The key that should be searched.
    void erase(const KeyT &key) {
        // TODO
    }

    /// Inserts a new entry into the tree.
    /// @param[in] key      The key that should be inserted.
    /// @param[in] value    The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        // TODO
    }
};

}  // namespace moderndbs

#endif
